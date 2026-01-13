import unittest
from bioinf.filogenia import distancia_levenshtein, upgma_simples

class TestFilogenia(unittest.TestCase):
    
    def test_levenshtein(self):
        # GATA vs ATA (Apaga G, Apaga T) -> Dist 2
        self.assertEqual(distancia_levenshtein("GATTA", "ATA"), 2)
        # Iguais -> 0
        self.assertEqual(distancia_levenshtein("ATG", "ATG"), 0)
        
    def test_upgma_basico(self):
        # A e B são iguais, devem juntar-se primeiro
        seqs = ["A", "A", "C"] 
        # Esperado: (('A', 'A'), 'C') ou variante
        arvore = upgma_simples(seqs)
        
        # Verificar se é um tuplo (estrutura de árvore)
        self.assertIsInstance(arvore, tuple)
        # Verificar se tem 2 ramos principais
        self.assertEqual(len(arvore), 2)

if __name__ == '__main__':
    unittest.main()