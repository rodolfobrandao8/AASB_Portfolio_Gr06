import unittest
from bioinf.blast import blast_simplificado

class TestBlast(unittest.TestCase):
    
    def test_blast_match_perfeito(self):
        # Encontrar "GATA" dentro de uma seq maior
        # A semente w=3 ("GAT" ou "ATA") deve disparar o hit
        db = "CCCGATACCC"
        query = "GATA"
        
        # Esperamos encontrar GATA (score 4*2 = 8) na posição 3
        q_res, db_res, score, start = blast_simplificado(query, db, w=3)
        
        self.assertEqual(db_res, "GATA")
        self.assertEqual(start, 3)
        self.assertEqual(score, 8)

    def test_blast_sem_match(self):
        # Se não houver semente comum, retorna vazio
        db = "AAAAAA"
        query = "TTTTTT"
        
        q_res, db_res, score, start = blast_simplificado(query, db, w=3)
        self.assertEqual(score, 0)
        
    def test_extensao(self):
        # A semente é "AAA" (w=3)
        # Mas a extensão deve apanhar o "G" antes e depois
        # Query:  GAAAG
        # DB:    CCGAAAGCC
        db = "CCGAAAGCC"
        query = "GAAAG"
        
        # Semente 'AAA' faz match. Extensão apanha os Gs.
        q_res, db_res, score, start = blast_simplificado(query, db, w=3)
        
        self.assertEqual(db_res, "GAAAG")
        self.assertEqual(score, 10) # 5 matches * 2

if __name__ == '__main__':
    unittest.main()