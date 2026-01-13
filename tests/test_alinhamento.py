import unittest
from bioinf.alinhamento import needleman_wunsch

class TestAlinhamento(unittest.TestCase):
    
    def test_alinhamento_simples(self):
        # Sequências iguais devem dar match perfeito
        # Score esperado: 4 matches * 2 = 8
        s1 = "ATGC"
        s2 = "ATGC"
        a1, a2, score = needleman_wunsch(s1, s2, match=2, mismatch=-1, gap=-1)
        
        self.assertEqual(a1, "ATGC")
        self.assertEqual(a2, "ATGC")
        self.assertEqual(score, 8)

    def test_alinhamento_com_gap(self):
        # Exemplo: ATC vs AC (deve inserir gap no T)
        # A(2) + T/Gap(-1) + C(2) = Score 3
        s1 = "ATC"
        s2 = "AC"
        a1, a2, score = needleman_wunsch(s1, s2, match=2, mismatch=-1, gap=-1)
        
        self.assertEqual(a1, "ATC")
        self.assertEqual(a2, "A-C")
        self.assertEqual(score, 3)

    def test_alinhamento_mismatch(self):
        # Exemplo: A vs T (Mismatch) -> Score -1
        s1 = "A"
        s2 = "T"
        a1, a2, score = needleman_wunsch(s1, s2, match=2, mismatch=-1, gap=-1)
        
        self.assertEqual(score, -1) 
        # Alinhamento pode ser A/T (mismatch) ou -A/T- (gaps), depende da preferencia
        # Com match=2, mismatch=-1, gap=-1, o mismatch (-1) é melhor que 2 gaps (-2)

    def test_smith_waterman_local(self):
        # Exemplo clássico: encontrar "TAG" dentro de "GGTAGCC"
        # Seq1: GGTAGCC
        # Seq2: TAG
        # O SW deve ignorar o início "GG" e o fim "CC" e alinhar só o "TAG"
        from bioinf.alinhamento import smith_waterman
        
        s1 = "GGTAGCC"
        s2 = "TAG"
        a1, a2, score = smith_waterman(s1, s2, match=2, mismatch=-1, gap=-1)
        
        self.assertEqual(a1, "TAG")
        self.assertEqual(a2, "TAG")
        self.assertEqual(score, 6) # 3 matches * 2

    def test_smith_waterman_sem_semelhanca(self):
        # Se não houver nada parecido, o score é 0 e strings vazias
        from bioinf.alinhamento import smith_waterman
        s1 = "AAAA"
        s2 = "TTTT"
        a1, a2, score = smith_waterman(s1, s2, match=2, mismatch=-5, gap=-5)
        
        self.assertEqual(score, 0)
        self.assertEqual(a1, "")
        self.assertEqual(a2, "")

    def test_alinhamento_multiplo(self):
        # Testar com 3 sequências que partilham um padrão "ATGC"
        # Seq1: ATGC
        # Seq2: AT-C (Falta o G)
        # Seq3: A-GC (Falta o T)
        # O consenso deve recuperar o máximo de informação possível
        from bioinf.alinhamento import alinhamento_multiplo
        
        seqs = ["ATGC", "ATC", "AGC"]
        consenso = alinhamento_multiplo(seqs)
        
        # O consenso ideal seria algo próximo de ATGC (comprimento 4)
        # Dependendo da ordem de fusão, pode variar, mas deve ter tamanho >= 3
        self.assertTrue(len(consenso) >= 3)
        # Verifica se as letras principais estão presentes
        self.assertIn("A", consenso)
        self.assertIn("C", consenso)

if __name__ == '__main__':
    unittest.main()