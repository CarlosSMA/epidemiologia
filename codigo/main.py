import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import os

def main():
    # --- 1. Parâmetros Epidemiológicos ---
    # População total fixa (equivalente ao total de células do grid 100x100)
    TAMANHO_POPULACAO = 10000

    # Beta: Taxa de transmissão (probabilidade de contágio)
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
    t = np.linspace(0, DIAS_TOTAIS, DIAS_TOTAIS)

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
    resultado = odeint(modelo_seir, y0, t, args=(TAMANHO_POPULACAO, BETA, SIGMA, GAMMA))
    S, E, I, R = resultado.T

    # --- 5. Geração do Relatório ---
    print("Gerando gráfico SEIR...")

    plt.figure(figsize=(10, 6))
    plt.plot(t, S, label='Suscetível (Saudável)', color='#31a354', linewidth=2.5)
    plt.plot(t, E, label='Exposto (Latente)', color='#fec44f', linewidth=2)
    plt.plot(t, I, label='Infectado (Transmitindo)', color='#de2d26', linewidth=3)
    plt.plot(t, R, label='Removido (Perda Foliar)', color='#969696', linewidth=2)

    plt.title('Dinâmica SEIR - Ferrugem Asiática (Abordagem Analítica)')
    plt.xlabel('Dias da Safra')
    plt.ylabel('Número de Plantas (População)')
    plt.legend()
    plt.grid(True)

    caminho_arquivo = 'relatorio/resultado_seir.png'
    plt.tight_layout()
    plt.savefig(caminho_arquivo, dpi=150)
    print(f"Gráfico salvo em: {caminho_arquivo}")

    plt.show()

if __name__ == "__main__":
    main()
