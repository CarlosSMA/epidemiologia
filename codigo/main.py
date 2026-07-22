import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import os

def main():
    # --- 1. Parâmetros Epidemiológicos ---
    TAMANHO_POPULACAO = 10_000
    BETA = 0.40

    # Sigma (Taxa de Latência): 1 / dias no estado Exposto (7 dias)
    SIGMA = 1.0 / 7.0

    # Gamma (Taxa de Remoção): 1 / dias no estado Infectado (14 dias)
    GAMMA = 1.0 / 14.0

    # Duração total da simulação
    DIAS_TOTAIS = 150

    # --- 2. Condições Iniciais ---
    I0 = 1                # Foco inicial da doença ("paciente zero")
    E0 = 0                # Nenhuma planta exposta no dia 0
    R0 = 0                # Nenhuma planta removida no dia 0
    S0 = TAMANHO_POPULACAO - I0 - E0 - R0 # Restante da lavoura começa sadia

    # Vetor de tempo (de 0 a 150 dias, com 150 pontos)
    tempo = np.linspace(0, DIAS_TOTAIS, DIAS_TOTAIS)

    # --- 3. Definição do Sistema de Equações Diferenciais ---
    def modelo_seir(y, t, N, beta, sigma, gamma):
        """
        Define as equações diferenciais que regem a transição entre os estados.
        """
        S, E, I, R = y

        # Taxas de variação diária para cada compartimento
        dSdt = -beta * S * I / N             # Suscetíveis diminuem à medida que se infectam
        dEdt = beta * S * I / N - sigma * E  # Expostos recebem contagiados e perdem os que viram Infectados
        dIdt = sigma * E - gamma * I         # Infectados recebem os Expostos e perdem os Removidos
        dRdt = gamma * I                     # Removidos acumulam os que deixam de ser Infectados

        return [dSdt, dEdt, dIdt, dRdt]

    # --- 4. Execução da Simulação ---
    print("Calculando a evolução da epidemia (EDOs)...")
    y0 = [S0, E0, I0, R0]
    # O método odeint resolve o sistema matemático para todos os dias de uma só vez
    resultado = odeint(modelo_seir, y0, tempo, args=(TAMANHO_POPULACAO, BETA, SIGMA, GAMMA))
    S, E, I, R = resultado.T

    criacao_relatorio(tempo, S, E, I, R)

def criacao_relatorio(tempo, S, E, I, R):
    print("Gerando gráfico SEIR...")
    if not os.path.exists('relatorio_ferrugem'):
        os.makedirs('relatorio_ferrugem')

    plt.figure(figsize=(10, 6))
    plt.plot(tempo, S, label='Suscetível (Saudável)', color='#31a354', linewidth=2.5)
    plt.plot(tempo, E, label='Exposto (Latente)', color='#fec44f', linewidth=2)
    plt.plot(tempo, I, label='Infectado (Transmitindo)', color='#de2d26', linewidth=3)
    plt.plot(tempo, R, label='Removido (Perda Foliar)', color='#969696', linewidth=2)

    plt.title('Dinâmica SEIR - Ferrugem Asiática (Abordagem Analítica)')
    plt.xlabel('Dias da Safra')
    plt.ylabel('Número de Plantas (População)')
    plt.legend()
    plt.grid(True)

    caminho_arquivo = 'relatorio_ferrugem/grafico_seir_edo.png'
    plt.tight_layout()
    plt.savefig(caminho_arquivo, dpi=150)
    print(f"Gráfico salvo em: {caminho_arquivo}")

    plt.show()

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
