import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from bioinf import sequencias, alinhamento, blast, filogenia, motifs

print("=== SEQUÊNCIAS ===")
dna = "ATGCGTAC"
rna = sequencias.transcricao(dna)
comp = sequencias.complemento(dna)
rev = sequencias.reverso(dna)
comp_inv = sequencias.complemento_inverso(dna)

print("DNA:", dna)
print("RNA:", rna)
print("Complemento:", comp)
print("Reverso:", rev)
print("Complemento reverso:", comp_inv)
print("Validar DNA:", sequencias.validar_dna(dna))
print("Validar RNA:", sequencias.validar_rna(rna))
print("Validar proteína 'ACDEFGHIKLMNPQRSTVWY':", sequencias.validar_proteina("ACDEFGHIKLMNPQRSTVWY"))

print("\n=== ALINHAMENTOS ===")
seq1 = "ACGTACGT"
seq2 = "ACGTCGTA"
a1, a2, score = alinhamento.needleman_wunsch(seq1, seq2)
print("NW Global Alinhamento:")
print(a1)
print(a2)
print("Score:", score)

a1_sw, a2_sw, score_sw = alinhamento.smith_waterman(seq1, seq2)
print("SW Local Alinhamento:")
print(a1_sw)
print(a2_sw)
print("Score:", score_sw)

seqs_multi = ["ACGT", "ACGG", "ACGA"]
alinhamento_multi, consenso = alinhamento.alinhamento_multiplo(seqs_multi)
print("Alinhamento múltiplo:", alinhamento_multi)
print("Consenso:", consenso)

print("\n=== BLAST SIMPLIFICADO ===")
query = "ACGTACGT"
target = "TTACGTACGTTT"
sub_q, sub_t, score_blast, pos = blast.blast_simplificado(query, target)
print("Query:", query)
print("Target:", target)
print("Melhor HSP Query:", sub_q)
print("Melhor HSP Target:", sub_t)
print("Score:", score_blast, "Posição:", pos)

print("\n=== FILOGENIA ===")
seqs_filo = ["ACGT", "ACGG", "ACGA", "TCGA"]
dist = filogenia.matriz_distancias(seqs_filo)
print("Matriz de distâncias:", dist)
arvore = filogenia.upgma(seqs_filo)
print("Árvore UPGMA simplificada:", arvore)

print("\n=== MOTIFS ===")
seq_motif = "GAATTCGAATTC"
padrao = "G^AATTC"
fragmentos, cortes = motifs.fragmentar_dna(seq_motif, padrao)
print("Sequência:", seq_motif)
print("Fragmentos:", fragmentos)
print("Cortes:", cortes)

lista_seqs = ["ACGT", "ACGG", "ACGA"]
pwm = motifs.criar_pwm(lista_seqs)
print("PWM:", pwm)

subseq_pos, subseq, prob = motifs.subsequencia_mais_provavel(pwm, "ACGTACTGAC")
print("Subsequência mais provável:", subseq, "Posição:", subseq_pos, "Probabilidade:", prob)

pssm = motifs.pwm_para_pssm(pwm)
score_pssm = motifs.score_seq_pssm(pssm, "ACGT")
print("PSSM:", pssm)
print("Score de sequência usando PSSM:", score_pssm)
