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

dna = "ATGCGT"
print(sequencias.validar_dna(dna))       # True
print(sequencias.transcricao(dna))       # AUGCGU
print(sequencias.complemento_inverso(dna)) # ACGCAT


```

## Alinhamento.py

```bash
from bioinf import alinhamento

seq1 = "ACGT"
seq2 = "AGT"
a1, a2, score = alinhamento.needleman_wunsch(seq1, seq2)
print(a1, a2, score)


```

## Blast.py

```bash
from bioinf import blast

query = "ATGC"
target = "TATGCATG"
sub_q, sub_t, score, pos = blast.blast_simplificado(query, target)
print(sub_q, sub_t, score, pos)


```
## Filogenia.py

```bash
from bioinf import filogenia

seqs = ["ATG", "ATC", "AGC"]
tree = filogenia.upgma(seqs)
print(tree)

```
## Motifs.py

```bash
from bioinf import motifs

seq = "ATGCGATG"
pattern = "A-x-G"
positions = motifs.procurar_motifs(seq, pattern)
print(positions)



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
