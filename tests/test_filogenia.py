import unittest
from bioinf.filogenia import distancia_levenshtein, matriz_distancias, upgma

class TestDistanciaLevenshtein(unittest.TestCase):
    def test_sequencias_identicas(self):
        """Distância entre sequências idênticas deve ser 0"""
        self.assertEqual(distancia_levenshtein("ATGC", "ATGC"), 0)

    def test_sequencias_completamente_diferentes(self):
        """Distância entre sequências totalmente diferentes"""
        self.assertEqual(distancia_levenshtein("AAAA", "TTTT"), 4)

    def test_sequencias_vazias(self):
        """Sequências vazias devem retornar distância correta"""
        self.assertEqual(distancia_levenshtein("", ""), 0)
        self.assertEqual(distancia_levenshtein("A", ""), 1)
        self.assertEqual(distancia_levenshtein("", "G"), 1)

    def test_sequencias_tamanho_1(self):
        """Sequências de tamanho 1"""
        self.assertEqual(distancia_levenshtein("A", "G"), 1)
        self.assertEqual(distancia_levenshtein("C", "C"), 0)


class TestMatrizDistancias(unittest.TestCase):
    def test_matriz_distancias_basica(self):
        """Matriz de distâncias par-a-par básica"""
        seqs = ["A", "G", "C"]
        m = matriz_distancias(seqs)
        self.assertEqual(m[("A", "G")], 1)
        self.assertEqual(m[("G", "A")], 1)  # simétrica
        self.assertEqual(m[("A", "C")], 1)
        self.assertEqual(m[("G", "C")], 1)

    def test_matriz_com_sequencia_vazia(self):
        """Matriz de distâncias com sequências vazias"""
        seqs = ["A", ""]
        m = matriz_distancias(seqs)
        self.assertEqual(m[("A", "")], 1)
        self.assertEqual(m[("", "A")], 1)


class TestUPGMA(unittest.TestCase):
    def test_upgma_basico(self):
        """UPGMA cria tupla aninhada para 3 sequências"""
        seqs = ["A", "G", "C"]
        arvore = upgma(seqs)
        # Deve ser uma tupla aninhada contendo as três sequências
        self.assertIn("A", str(arvore))
        self.assertIn("G", str(arvore))
        self.assertIn("C", str(arvore))

    def test_upgma_com_sequencia_vazia(self):
        """UPGMA com sequência vazia"""
        seqs = ["", "A"]
        arvore = upgma(seqs)
        self.assertIn("", str(arvore))
        self.assertIn("A", str(arvore))

    def test_upgma_tamanho_1(self):
        """UPGMA com apenas uma sequência"""
        seqs = ["A"]
        arvore = upgma(seqs)
        self.assertEqual(arvore, "A")


if __name__ == "__main__":
    unittest.main()
