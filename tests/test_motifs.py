import unittest
from bioinf.motifs import prosite_para_regex, procurar_motifs, criar_pwm, subsequencia_mais_provavel

class TestMotifs(unittest.TestCase):
    
    def test_prosite_conversion(self):
        # Testar a regra: x -> . e {ED} -> [^ED] e (4) -> {4}
        # Exemplo da tua imagem do caderno
        prosite = "[AC]-x-V-x(4)-{ED}"
        regex = prosite_para_regex(prosite)
        self.assertEqual(regex, "[AC].V.{4}[^ED]")
        
    def test_procurar_motif(self):
        # O motif C-A-T deve ser encontrado em GGCATGG
        # Index esperado: 2
        seq = "GGCATGG"
        padrao = "C-A-T"
        pos = procurar_motifs(seq, padrao)
        self.assertEqual(pos, [2])
        
    def test_pwm_laplace(self):
        # Testar cálculo simples com pseudocount=1
        # Seqs: A, A. (2 sequências)
        # Pos 1: 2 As. Laplace: (2+1)/(2+4) = 3/6 = 0.5
        seqs = ["A", "A"]
        pwm = criar_pwm(seqs)
        
        # A probabilidade de A na posição 0 deve ser 0.5
        self.assertAlmostEqual(pwm[0]['A'], 0.5)
        # A probabilidade de C (que não existe) deve ser (0+1)/6 = 0.166...
        self.assertAlmostEqual(pwm[0]['C'], 1/6)

    def test_subsequencia_provavel(self):
        # Se a PWM prefere A, deve escolher AAA em vez de TTT
        seqs = ["AAA", "AAA"]
        pwm = criar_pwm(seqs)
        
        alvo = "TTTAAATTT"
        pos, sub, prob = subsequencia_mais_provavel(pwm, alvo)
        
        self.assertEqual(sub, "AAA")
        self.assertEqual(pos, 3) # Começa no índice 3

if __name__ == '__main__':
    unittest.main()