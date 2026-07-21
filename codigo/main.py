from resultados import Resultado
import numpy as np
import random
from matplotlib.colors import ListedColormap
from tqdm import tqdm

# --- Configurações da Simulação ---
tamanho = 100             # Tamanho da lavoura (NxN)
numero_dias = 150    # Duração da simulação em dias

# Parâmetros Epidemiológicos (Ajustados para Ferrugem)
# Probabilidade de uma planta doente infectar uma vizinha sadia em um dia.
BETA = 0.40
DAYS_LATENT = 7     # Dias entre o contágio e a manifestação dos sintomas (E -> I)
DAYS_INFECTIOUS = 14 # Dias que a planta libera esporos antes de morrer/ser removida (I -> R)

# Configurações de Relatório
IMAGE_INTERVAL = 15  # Grava imagem do grid a cada 15 dias

# --- Inicialização ---
def initialize_grids(n_size):
    # Grade de estados principais (0, 1, 2, 3, 4)
    state_grid = np.zeros((n_size, n_size), dtype=int)
    # Grades temporizadoras para contar os dias em cada estado
    latent_timers = np.zeros((n_size, n_size), dtype=int)
    infectious_timers = np.zeros((n_size, n_size), dtype=int)
    return state_grid, latent_timers, infectious_timers

# --- Função de Plantio Aleatório ---
def run_planting_phase(state_grid, n_size):
    # Requisitos: Garantir que o total de soja (S+E+I+R) esteja entre 30% e 80% do campo (10000 células)
    MIN_SOYBEAN_COUNT = 3000
    MAX_SOYBEAN_COUNT = 8000

    # 1. Conta o estado atual
    soybean_mask = (state_grid >= 1)
    current_soybean_count = np.sum(soybean_mask)

    # Define um alvo de população aleatório para hoje
    target_population = random.randint(MIN_SOYBEAN_COUNT, MAX_SOYBEAN_COUNT)

    # 2. Se a população atual está abaixo do alvo, planta novas S (estado 1)
    if current_soybean_count < target_population:
        # Define quantas plantas precisam ser adicionadas
        needed_plants = target_population - current_soybean_count

        # Encontra os índices de todas as células vazias (estado 0)
        empty_indices = np.argwhere(state_grid == 0)

        if len(empty_indices) > 0:
            # Seleciona aleatoriamente onde plantar
            plants_to_add = min(needed_plants, len(empty_indices))
            random_indices = np.random.choice(len(empty_indices), plants_to_add, replace=False)
            chosen_locations = empty_indices[random_indices]

            # Atualiza o estado
            state_grid[chosen_locations[:, 0], chosen_locations[:, 1]] = 1

# --- Lógica de Contágio (Regra de Vizinhança) ---
def calculate_infection_probability(grid, r, c, n_size, beta_rate):
    # Procura vizinhos infectados (estado 3) na Vizinhança de Moore (8 vizinhos)

    # Tratamento de bordas (não-toroidal: as plantas na borda têm menos vizinhos)
    min_r, max_r = max(0, r-1), min(n_size-1, r+1)
    min_c, max_c = max(0, c-1), min(n_size-1, c+1)

    # Conta quantos vizinhos estão no estado 3 (Infectados)
    infected_neighbors = 0
    for ir in range(min_r, max_r + 1):
        for ic in range(min_c, max_c + 1):
            if (ir == r and ic == c): continue # Não conta a própria célula
            if grid[ir, ic] == 3:
                infected_neighbors += 1

    if infected_neighbors == 0:
        return 0.0

    # Probabilidade total de contágio se ao menos um vizinho estiver doente
    # Probabilidade = 1 - (1 - beta) ^ numero_de_vizinhos_doentes
    return 1 - (1 - beta_rate) ** infected_neighbors

# --- Execução da Simulação ---
state_grid, latent_timers, infectious_timers = initialize_grids(tamanho)

# Adiciona o foco inicial da doença ("Paciente Zero") no centro
state_grid[tamanho//2, tamanho//2] = 3

# Histórico da população para gráficos
historico_populacao = {'S': [], 'E': [], 'I': [], 'R': [], 'Empty': []}

print(f"Iniciando simulação de Ferrugem Asiática (SEIR) por {numero_dias} dias...")

# Loop principal de tempo (cada tick = 1 dia)
for day in tqdm(range(1, numero_dias + 1)):

    # 1. Fase de Plantio Dinâmico
    run_planting_phase(state_grid, tamanho)

    # 2. Criar uma cópia da grade para gravar o próximo estado (Double Buffering)
    next_state = np.copy(state_grid)
    next_latent = np.copy(latent_timers)
    next_infectious = np.copy(infectious_timers)

    # 3. Aplicar as regras de transição para cada célula
    for r in range(tamanho):
        for c in range(tamanho):
            current_cell_state = state_grid[r, c]

            # --- Regras para plantas Suscetíveis (1) ---
            if current_cell_state == 1:
                prob = calculate_infection_probability(state_grid, r, c, tamanho, BETA)
                if prob > 0 and random.random() < prob:
                    next_state[r, c] = 2 # Transição S -> Exposto
                    next_latent[r, c] = 1 # Inicia o timer latente

            # --- Regras para plantas Expostas (Latentes) (2) ---
            elif current_cell_state == 2:
                next_latent[r, c] += 1
                # Fim do período de latência: Transição E -> Infectado
                if next_latent[r, c] >= DAYS_LATENT:
                    next_state[r, c] = 3
                    next_latent[r, c] = 0 # Zera o timer latente
                    next_infectious[r, c] = 1 # Inicia o timer de infecção

            # --- Regras para plantas Infectadas (Ativas) (3) ---
            elif current_cell_state == 3:
                next_infectious[r, c] += 1
                # Fim do período de infecção: Transição I -> Removido
                if next_infectious[r, c] >= DAYS_INFECTIOUS:
                    next_state[r, c] = 4
                    next_infectious[r, c] = 0 # Zera o timer de infecção

            # --- Regras para plantas Removidas (4) ---
            # Neste modelo simplificado, as plantas Removidas não reaparecem como solo vazio (0).
            # Para simular perda total da folha, você poderia adicionar: I -> 0.
            # Manteremos I->R para traçar o total acumulado da doença no gráfico SEIR.
            elif current_cell_state == 4:
                pass

    # Atualiza as grades principais para o próximo dia
    state_grid, latent_timers, infectious_timers = next_state, next_latent, next_infectious

    # 4. Gravar dados e instantâneos
    counts = np.bincount(state_grid.flatten(), minlength=5)
    historico_populacao['Empty'].append(counts[0])
    historico_populacao['S'].append(counts[1])
    historico_populacao['E'].append(counts[2])
    historico_populacao['I'].append(counts[3])
    historico_populacao['R'].append(counts[4])

    # Grava imagem do grid em intervalos
    if day % IMAGE_INTERVAL == 0 or day == 1:
        grid_snapshots.append(np.copy(state_grid))
        day_labels.append(day)

resultado = Resultado(historico_populacao, numero_dias, tamanho)

print("\nSimulação concluída. Gerando relatório final...")
