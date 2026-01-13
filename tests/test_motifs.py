import unittest
import math
from bioinf import motifs


class TestPrositeParaRegex(unittest.TestCase):
    def test_regex_basico(self):
        padrao = "A-x-{C}-G(2,3)"
        regex = motifs.prosite_para_regex(padrao)
        self.assertIn("A", regex)
        self.assertIn("[^C]", regex)
        self.assertIn("{2,3}", regex)

    def test_ancoras(self):
        padrao = "<A-CG>"
        regex = motifs.prosite_para_regex(padrao)
        self.assertTrue(regex.startswith("^"))
        self.assertTrue(regex.endswith("$"))


class TestProcurarMotifs(unittest.TestCase):
    def test_motifs_simples(self):
        seq = "ATGCGATG"
        padrao = "ATG"
        posicoes = motifs.procurar_motifs(seq, padrao)
        self.assertEqual(posicoes, [0, 5])

    def test_sem_ocorrencias(self):
        seq = "AAAA"
        padrao = "TTT"
        posicoes = motifs.procurar_motifs(seq, padrao)
        self.assertEqual(posicoes, [])


class TestEnzimaFragmentacao(unittest.TestCase):
    def test_enzima_para_regex(self):
        seq, corte = motifs.enzima_para_regex("G^AATTC")
        self.assertEqual(seq, "GAATTC")
        self.assertEqual(corte, 1)

    def test_fragmentar_dna_basico(self):
        seq = "GAATTCCGAATT"
        frag, cortes = motifs.fragmentar_dna(seq, "G^AATTC")
        self.assertEqual(frag, ["G", "CGAATT", ""])
        self.assertEqual(cortes, [1, 7])


class TestPWM(unittest.TestCase):
    def test_criar_pwm_simples(self):
        seqs = ["ACG", "ACG"]
        pwm_result = motifs.criar_pwm(seqs, pseudocount=0)
        self.assertAlmostEqual(pwm_result[0]["A"], 1.0)
        self.assertAlmostEqual(pwm_result[1]["C"], 1.0)
        self.assertAlmostEqual(pwm_result[2]["G"], 1.0)

    def test_probabilidade_seq_pwm(self):
        seqs = ["ACG"]
        pwm_result = motifs.criar_pwm(seqs)
        prob = motifs.probabilidade_seq_pwm(pwm_result, "ACG")
        self.assertGreater(prob, 0)

    def test_subsequencia_mais_provavel(self):
        seqs = ["ACG", "ACG"]
        pwm_result = motifs.criar_pwm(seqs)
        pos, sub, prob = motifs.subsequencia_mais_provavel(pwm_result, "TTACGTT")
        self.assertEqual(sub, "ACG")
        self.assertEqual(pos, 2)

    def test_pwm_para_pssm_e_score(self):
        seqs = ["ACG"]
        pwm_result = motifs.criar_pwm(seqs)
        pssm = motifs.pwm_para_pssm(pwm_result)
        score = motifs.score_seq_pssm(pssm, "ACG")
        self.assertGreater(score, 0)
        score_bad = motifs.score_seq_pssm(pssm, "TTT")
        self.assertLess(score_bad, score)


class TestCasosLimite(unittest.TestCase):
    def test_lista_vazia_pwm(self):
        pwm_result = motifs.criar_pwm([])
        self.assertEqual(pwm_result, [])

    def test_sequencia_vazia_probabilidade(self):
        pwm_result = motifs.criar_pwm(["A"])
        prob = motifs.probabilidade_seq_pwm(pwm_result, "")
        self.assertEqual(prob, 1.0)


if __name__ == "__main__":
    unittest.main()
