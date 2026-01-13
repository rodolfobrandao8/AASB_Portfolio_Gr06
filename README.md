# Portefólio de Algoritmos para análise de sequências biológicas (AASB 2025/2026)

Biblioteca Python que implementa os principais algoritmos abordados na UC.
O projeto foi desenvolvido com foco em correção algorítmica, qualidade de código, documentação completa e testes unitários, sendo facilmente importável e reutilizável por terceiros.

---

## Conteúdo do Projeto

Este portefólio inclui implementações dos seguintes tópicos:

### 1. Sequências Biológicas
- Validação de sequências (DNA, RNA e proteínas)
- Complemento e reverso-complemento
- Transcrição DNA → RNA

### 2. Alinhamento de Sequências
- Matrizes de substituição (PAM / BLOSUM)
- Needleman–Wunsch (alinhamento global)
- Smith–Waterman (alinhamento local)
- Alinhamento múltiplo progressivo
- Reconstrução dos alinhamentos e cálculo de score

### 3. Motifs e Padrões
- Procura de padrões fixos com ambiguidades
- Conversão PROSITE → expressões regulares
- Enzimas de restrição → regex + posições de corte
- PWM e PSSM:
  - Construção
  - Probabilidade de uma sequência
  - Subsequência mais provável

### 4. BLAST Simplificado
- Query map
- Identificação de hits
- Extensão dos alinhamentos
- Melhor alinhamento local

### 5. Análise Filogenética
- Construção de matriz de distâncias
- UPGMA:
  - Clustering
  - Construção da árvore filogenética


# Uso dos códigos

## Sequencias.py

```bash
from bioinf import sequencias

# Validação
sequencias.validar_dna("ACGT")       # True
sequencias.validar_rna("ACGU")       # True
sequencias.validar_proteina("ACDE")  # True

# Transcrição e complementos
rna = sequencias.transcricao("ACGT")            # ACGU
comp = sequencias.complemento("ACGT")          # TGCA
rev = sequencias.reverso("ACGT")               # TGCA
comp_inv = sequencias.complemento_inverso("ACGT") # ACGT



```

## Alinhamento.py

```bash
from bioinf import alinhamento

seq1 = "ACGT"
seq2 = "AGT"

# Alinhamento global (Needleman-Wunsch)
a1, a2, score = alinhamento.needleman_wunsch(seq1, seq2)
print("Global:", a1, a2, score)

# Alinhamento local (Smith-Waterman)
a1, a2, score = alinhamento.smith_waterman(seq1, seq2)
print("Local:", a1, a2, score)

# Alinhamento múltiplo
seqs = ["ACGT", "AGT", "ACG"]
alns, consenso = alinhamento.alinhamento_multiplo(seqs)
print("Múltiplo:", alns, "Consenso:", consenso)

# Dot plot
dot = alinhamento.dot_plot(seq1, seq2)
print("Dot plot:", dot)

# Consenso de múltiplas sequências
cons = alinhamento.consenso_multiplas(["ACG", "AGG", "ACG"])
print("Consenso:", cons)



```

## Blast.py

```bash
from bioinf import blast

query = "ACGTAC"
seq_alvo = "TTACGTACGG"

sub_q, sub_t, score, start = blast.blast_simplificado(query, seq_alvo)
print(sub_q, sub_t, score, start)



```
## Filogenia.py

```bash
from bioinf import filogenia

seqs = ["ACGT", "AGT", "ACG"]

# Distância de Levenshtein
d = filogenia.distancia_levenshtein("ACGT", "AGT")
print("Distância:", d)

# Matriz de distâncias
matriz = filogenia.matriz_distancias(seqs)
print("Matriz de distâncias:", matriz)

# Árvore filogenética simplificada (UPGMA)
arvore = filogenia.upgma(seqs)
print("Árvore:", arvore)


```
## Motifs.py

```bash
from bioinf import motifs

# PROSITE → regex
regex = motifs.prosite_para_regex("A-x-{C}-G(2,3)")

# Procurar motifs
posicoes = motifs.procurar_motifs("ATGCGATG", "ATG")
print("Posições:", posicoes)

# Fragmentação de DNA
fragmentos, cortes = motifs.fragmentar_dna("GAATTCCGAATT", "G^AATTC")
print("Fragmentos:", fragmentos, "Cortes:", cortes)

# PWM/PSSM
seqs = ["ACG", "ACG"]
pwm = motifs.criar_pwm(seqs)
pos, sub, prob = motifs.subsequencia_mais_provavel(pwm, "TTACGTT")
pssm = motifs.pwm_para_pssm(pwm)
score = motifs.score_seq_pssm(pssm, "ACG")




```
# Testes Unitários
Os testes foram desenvolvidos com unittest / pytest e cobrem:
Casos normais
Casos limite (sequências vazias, tamanho 1)
Exceções
Comparações com exemplos da literatura

## Executar testes

```bash
pytest
```

## Executar testes com cobertura
``` bash
pytest --cov=bioinf --cov-report=term-missing
```

# Documentação (Sphinx)
A documentação inclui:
Página inicial
Guia de instalação
Tutorial de utilização
Referência completa da API (docstrings)

# Qualidade do Código (Radon)
O código foi desenvolvido de forma a garantir:
Complexidade ciclomática aceitável
Funções pequenas e legíveis
Conformidade com PEP 8

## Avaliar complexidade:

```bash
radon cc bioinf/ -a -s
radon mi bioinf/ -s
```

# Autores
Grupo 06
UC: Algoritmos e Análise de Sistemas Biológicos
Ano letivo: 2025/2026
