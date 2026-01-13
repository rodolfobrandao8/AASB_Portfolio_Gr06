import unittest
from bioinf.alinhamento import (
    needleman_wunsch,
    smith_waterman,
    dot_plot,
    consenso_multiplas,
    alinhamento_multiplo,
    BLOSUM62
)

class TestNeedlemanWunsch(unittest.TestCase):
    def test_sequencias_identicas(self):
        a1, a2, score = needleman_wunsch("ATGC", "ATGC", BLOSUM62)
        self.assertEqual(a1, "ATGC")
        self.assertEqual(a2, "ATGC")
        self.assertGreater(score, 0)

    def test_sequencia_vazia(self):
        a1, a2, score = needleman_wunsch("", "", BLOSUM62)
        self.assertEqual(a1, "")
        self.assertEqual(a2, "")
        self.assertEqual(score, 0)

    def test_sequencias_tamanho_1(self):
        a1, a2, score = needleman_wunsch("A", "T", BLOSUM62)
        self.assertEqual(len(a1), 1)
        self.assertEqual(len(a2), 1)

class TestSmithWaterman(unittest.TestCase):
    def test_alinhamento_local_simples(self):
        a1, a2, score = smith_waterman("ATGC", "TGC", BLOSUM62)
        self.assertIn("TGC", a1+a2)
        self.assertGreaterEqual(score, 0)

    def test_sequencia_vazia(self):
        a1, a2, score = smith_waterman("", "", BLOSUM62)
        self.assertEqual(a1, "")
        self.assertEqual(a2, "")
        self.assertEqual(score, 0)

class TestDotPlot(unittest.TestCase):
    def test_dotplot_basico(self):
        matriz = dot_plot("AT", "AG")
        self.assertEqual(matriz, [[1,0],[0,0]])

class TestConsensoMultiplo(unittest.TestCase):
    def test_consenso_simples(self):
        alin = ["ATGC", "ATGA", "ATGT"]
        c = consenso_multiplas(alin)
        self.assertEqual(c, "ATGA")

class TestAlinhamentoMultiplo(unittest.TestCase):
    def test_alinhamento_multiplo_basico(self):
        seqs = ["ATGC", "ATGA", "ATGT"]
        alin, cons = alinhamento_multiplo(seqs, BLOSUM62)
        self.assertEqual(cons, "ATGA")
        self.assertTrue(all(len(a) == len(alin[0]) for a in alin))

    def test_sequencia_unica(self):
        seqs = ["A"]
        alin, cons = alinhamento_multiplo(seqs, BLOSUM62)
        self.assertEqual(alin, [["A"]])
        self.assertEqual(cons, "A")

    def test_sequencias_vazias(self):
        seqs = ["", "", ""]
        alin, cons = alinhamento_multiplo(seqs, BLOSUM62)
        self.assertTrue(all(a == "" or a == [""] for a in alin))
        self.assertEqual(cons, "")

if __name__ == "__main__":
    unittest.main()
