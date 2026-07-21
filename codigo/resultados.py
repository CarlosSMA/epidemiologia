import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os

class Resultado:
    def __init__(self, numero_dias, historico_populacao, tamanho):
        self.numero_dias = numero_dias
        self.historico_populacao = historico_populacao
        self.tamanho = tamanho

    def criar_relatorio(self):
        if not os.path.exists('relatorio_ferrugem'):
            os.makedirs('relatorio_ferrugem')

        # 1. Gráficos SEIR (Linha do tempo global)
        days_range = range(1, self.numero_dias + 1)
        plt.figure(figsize=(12, 7))
        plt.plot(days_range, self.historico_populacao['S'], label='Suscetível (Saudável)', color='#31a354', linewidth=2.5) # Verde escuro
        plt.plot(days_range, self.historico_populacao['E'], label='Exposto (Latente)', color='#fec44f', linewidth=2)    # Amarelo/Laranja
        plt.plot(days_range, self.historico_populacao['I'], label='Infectado (Transmitindo)', color='#de2d26', linewidth=3) # Vermelho
        plt.plot(days_range, self.historico_populacao['R'], label='Removido (Perda Foliar)', color='#969696', linewidth=2)  # Cinza

        plt.title(f'Relatório Populacional SEIR - Proliferação da Ferrugem Asiática (Total: {self.N*self.N} cells)')
        plt.xlabel('Dias da Safra (Tick)')
        plt.ylabel('Contagem de Plantas')
        plt.legend()
        plt.grid(True)
        plt.savefig('relatorio_ferrugem/grafico_seir_total.png', dpi=150)
        plt.close()


        # 2. Relatório de Imagens (snapshots do grid)
        # Mapa de cores personalizado (0=Vazio, 1=S, 2=E, 3=I, 4=R)
        cmap = ListedColormap(['#eaeaea', '#31a354', '#fec44f', '#de2d26', '#969696'])

        grid_snapshots = []
        day_labels = []

        num_snapshots = len(grid_snapshots)
        fig, axes = plt.subplots(1, num_snapshots, figsize=(5 * num_snapshots, 5))
        plt.suptitle('Evolução Espacial da Ferrugem Asiática na Lavoura (Autômato Celular)', fontsize=16)

        # Trata o caso de apenas um snapshot (axes não é uma lista)
        if num_snapshots == 1: axes = [axes]

        for ax, snap, day_num in zip(axes, grid_snapshots, day_labels):
            ax.imshow(snap, cmap=cmap, origin='lower', interpolation='nearest')
            ax.set_title(f'Dia {day_num}')
            ax.set_xticks([]) # Remove marcações de eixo
            ax.set_yticks([])

        plt.tight_layout(rect=(0, 0.03, 1, 0.95))
        plt.savefig('relatorio_ferrugem/relatorio_imagens_grid.png', dpi=150)
        plt.close()

        # Mostra o relatório final
        relatorio_imagem = plt.imread('relatorio_ferrugem/relatorio_imagens_grid.png')
        plt.figure(figsize=(16, 8))
        plt.imshow(relatorio_imagem)
        plt.axis('off')
        plt.show()

        relatorio_grafico = plt.imread('relatorio_ferrugem/grafico_seir_total.png')
        plt.figure(figsize=(16, 9))
        plt.imshow(relatorio_grafico)
        plt.axis('off')
        plt.show()

        print("\nRelatório gerado com sucesso na pasta 'relatorio_ferrugem'.")
