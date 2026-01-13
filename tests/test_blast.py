import unittest
from bioinf.blast import (
    blast_simplificado,
    construir_mapa,
    encontrar_hits,
    estender_hit
)


class TestConstruirMapa(unittest.TestCase):
    def test_mapa_basico(self):
        query = "ATGCAT"
        mapa = construir_mapa(query, 3)
        self.assertIn("ATG", mapa)
        self.assertIn("TGC", mapa)
        self.assertEqual(mapa["ATG"], [0])
    
    def test_query_tamanho_menor_w(self):
        query = "AT"
        mapa = construir_mapa(query, 3)
        self.assertEqual(mapa, {})


class TestEncontrarHits(unittest.TestCase):
    def test_hit_simples(self):
        seq = "ATGCAT"
        mapa = {"ATG": [0]}
        hits = encontrar_hits(seq, mapa, 3)
        self.assertIn((0, 0), hits)

    def test_sem_hits(self):
        seq = "AAAAA"
        mapa = {"TTT": [0]}
        hits = encontrar_hits(seq, mapa, 3)
        self.assertEqual(hits, [])


class TestEstenderHit(unittest.TestCase):
    def test_hit_basico(self):
        query = "ATGC"
        seq = "ATGC"
        hit = (0, 0)
        score, q_start, t_start, length = estender_hit(query, seq, hit, 3, match=2, mismatch=-1)
        self.assertEqual(q_start, 0)
        self.assertEqual(t_start, 0)
        self.assertEqual(length, 4)
        self.assertGreater(score, 0)

    def test_hit_com_mismatch(self):
        query = "ATGC"
        seq = "ATGA"
        hit = (0, 0)
        score, _, _, length = estender_hit(query, seq, hit, 3, match=2, mismatch=-1)
        self.assertEqual(length, 4)


class TestBlastSimplificado(unittest.TestCase):
    def test_alinhamento_perfeito(self):
        query = "ATGC"
        seq = "ATGC"
        sub_q, sub_t, score, pos = blast_simplificado(query, seq, w=2)
        self.assertEqual(sub_q, query)
        self.assertEqual(sub_t, seq)
        self.assertEqual(pos, 0)
        self.assertGreater(score, 0)

    def test_subsequencia_no_alvo(self):
        query = "ATGC"
        seq = "TTATGCGG"
        sub_q, sub_t, score, pos = blast_simplificado(query, seq, w=2)
        self.assertEqual(sub_q, "ATGC")
        self.assertEqual(sub_t, "ATGC")
        self.assertGreater(score, 0)

    def test_sem_hits(self):
        query = "AAAA"
        seq = "TTTT"
        sub_q, sub_t, score, pos = blast_simplificado(query, seq, w=2)
        self.assertEqual(sub_q, "")
        self.assertEqual(sub_t, "")
        self.assertEqual(score, 0)
        self.assertEqual(pos, -1)

    def test_query_menor_w(self):
        query = "AT"
        seq = "ATGC"
        sub_q, sub_t, score, pos = blast_simplificado(query, seq, w=3)
        self.assertEqual(sub_q, "")
        self.assertEqual(sub_t, "")
        self.assertEqual(score, 0)
        self.assertEqual(pos, -1)


if __name__ == "__main__":
    unittest.main()
