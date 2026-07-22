# Modelagem SEIR para Análise da Proliferação da Ferrugem Asiática em Plantações de Soja

Este repositório contém o código-fonte referente ao estudo de modelagem matemática da proliferação da Ferrugem Asiática na sojicultura brasileira. O projeto implementa um modelo epidemiológico compartimental determinístico **SEIR** (Suscetíveis, Expostas, Infectadas e Removidas) baseado em Equações Diferenciais Ordinárias (EDOs).

## Dependências
O projeto é desenvolvido em Python e depende das seguintes bibliotecas:
- `numpy`: Manipulação numérica
- `scipy`: Resolução do sistema de equações (através de `odeint`)
- `matplotlib`: Geração de gráficos

## Instruções de uso

1. Crie e ative um ambiente virtual Python (recomendado):
```shell
python -m venv .venv
source .venv/bin/activate
```

2. Instale as dependências necessárias:
```shell
pip install -r requirements.txt
```

3. Execute o arquivo principal:
```shell
python3 ./codigo/main.py
```

## Resultados
A execução do algoritmo criará um diretório `relatorio_ferrugem` na raiz do projeto e salvará o gráfico `grafico_seir_edo.png`. Este gráfico demonstra as curvas compartimentais da lavoura no decorrer da safra.
