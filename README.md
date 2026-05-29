
mle_fit
=======

Este projeto implementa o método desenvolvido por [Clauset et al. 2007](http://arxiv.org/abs/0706.1062) para ajuste de distribuições em cauda longa. Ele foi elaborado como parte do trabalho de monografia [Complexidade e Emergência em Linguagens Naturais](https://drive.google.com/file/d/1VF6BgdaVsb-DG5UdC0VnKhsgc-dNbC1J/view?usp=sharing) da Universidade Federal de Viçosa.

## ⚠️ Status do Projeto: Estágio Inicial

Atualmente, o `mle_fit` é um projeto em **desenvolvimento ativo**. A instalação e o uso ainda são manuais e serão aprimorados em breve. A API é funcional, mas será refinada para se tornar mais intuitiva.

## Funcionalidades

- **Ajuste em lei de potência:** Rotinas Fotran otimizadas para ajustes em alta-velocidade.
- **Cálculo estatística KS:** Rotina Fortran para cálculo da estatística de Kolmogorov-Smirnov (KS) de leis de potência.
- **Goodness-of-fit:** Cálculo do p-valor para determinar a plausibilidade da hipótese da lei de potência.
- **Plotting integrado:** Funções para traçar e customizar plots da distribuição, PDF e CCDF, com `matplotlib`.
- **Escolha de mínimos:** Permite a escolha de mínimos não-globais da estatística KS.

## Instalação (Manual)

**Pré-requisitos:**
- Compilador Fortran (ex: `gfortran`) instalado.
- Python 3.x.

Para instalar,
1. **Clone o repositório:**
```zsh
git clone https://github.com/timotheosf/mle_fit.git
cd mle_fit
```
2. **Compile os módulos em Fortran:**
```zsh
gfortran scr/_kstest.f90 -o kstest.out && gfortran scr/mle_discrete.f90 -o mle_discrete.out 
```
3. **Adicione a pasta ao PYTHONPATH:** 
```zsh 
export PYTHONPATH="${PYTHONPATH}:$(pwd)
```
*Nota: Este comando é temporário. Para torná-lo permanente, adicione-o ao seu arquivo de configuração do shell (ex: `.bashrc`, `.zshrc`).*

## Como usar

Depois da instalação, importe o módulo em seu script Python:
```python
import mle
```

### Ajuste de leis de potência

Para o uso básico de ajuste, basta fazer:
```python
data = array([ 1 , 2 , 3 , ... ]) #> Os dados devem ser um array
result = mle.fit_from_data( data )
x_min , alpha , std_alpha , ks_stats = result[0]
```
É importante mencionar que `result` é um array que contém várias outras informações:
```python
r_1 = result[0] #> devolve os parâmetros do ajuste no mínimo global
r_2 = result[1] #> devolve a lista de tuplas (x_min , alpha , std_alpha , ks_stats) para cada mínimo local de KS
r_3 = result[2] #> devolve a lista dos valores de KS para cada candidato a x_min
```


### Cálculo da estatística KS

Caso queira-se medir a estatística KS de um determinado conjunto de dados com o limitante inferior (`x_min`):
```python
ks_stats = mle.ks_stats_from_data( data , x_min )
```

### p-Valor

Para o cálculo do p-valor:
```python
p = mle.p_value_stats( data , x_min , alpha , ks , eps=0.01 )
```
onde `eps` representa a precisão do cálculo.

### Plot gráfico

Os comando para plotar a ccdf e a pdf das distribuições têm a mesma sintaxe:
```python
import matplotlib.pyplot as plt
#> PDF
ax = plt.gca( )
mle.pdf_plot( data , x_min , alpha , ax=ax , var=r'$\omega$' , expoent=r'\beta' )
# var é o nome da variável (no caso, omega) que ira apararecer, e expoent, o nome do expoente (alpha, beta, gamma, ...)
#> CDF
ax = plt.gca( )
mle.cdf_plot( data , x_min , alpha , ax=ax , var=r'$\omega$' , expoent=r'\beta' )
```

A log-verossimilhança e a estatística KS também podem ser visualizadas:
```python
#> Log-verossimilhança:
ax = plt.gca( )
mle.likelihood_plot( data , x_min , alpha , ax=ax )
#> Estatística KS (por mínimo)
ax = plt.gca( )
mle.ks_stats_plot( ks_stats=result[2] , ax=ax )
```

### Mínimos da estatística KS

Um dos diferenciais do pacote é permitir a escolha do mínimo da estatística KS. Em alguns casos, interessa que o ajuste capture o maior conjunto de dados em que ainda valha uma lei de potência ($p>0.1$), o que nem sempre coincide com a escolha do mínimo global. Para resolver esta ambiguidade, a função `fit_from_data` devolve uma lista com as informações necessárias para isso:
```python
result = mle.fit_from_data( data )
# Para o penúltimo mínimo:
x_min , alpha , std_alpha , ks_stats = result[1][-1][:]
# Para o antepenúltimo mínimo:
x_min , alpha , std_alpha , ks_stats = result[1][-2][:]
# Para o N-ésimo mínimo:
x_min , alpha , std_alpha , ks_stats = result[1][N][:]
```
Para auxiliar nessa escolha, a `ks_stats_plot` mostra a localização de todos os mínimos locais:
```python
ax = plt.gca( )
mle.ks_stats_plot( ks_stats=result[2] , ks_mins=result[1] , ax=ax , MINS=True , get_minimun=-1  )
# Aqui, get_minimun determina a posição da linha vertical do mínimo
```
