def validar_dna(seq):
    """Valida se uma sequência contém apenas nucleótidos de DNA.

    Aceita as bases ``A, C, G, T`` (case-insensitive). Se a sequência for vazia
    (``""``/``None``/falsy), devolve ``False``.

    Args:
        seq (str): Sequência a validar.

    Returns:
        bool: ``True`` se a sequência for não vazia e contiver apenas A/C/G/T;
        caso contrário ``False``.

    Raises:
        AttributeError: Se ``seq`` não for string e não tiver o método ``upper()``.

    Examples:
        >>> validar_dna("ACGT")
        True
        >>> validar_dna("acgt")
        True
        >>> validar_dna("ACGU")
        False
        >>> validar_dna("")
        False
    """
    if not seq:
        return False
    seq = seq.upper()
    return set(seq).issubset({'A', 'C', 'G', 'T'})


def validar_rna(seq):
    """Valida se uma sequência contém apenas nucleótidos de RNA.

    Aceita as bases ``A, C, G, U`` (case-insensitive). Se a sequência for vazia
    (``""``/``None``/falsy), devolve ``False``.

    Args:
        seq (str): Sequência a validar.

    Returns:
        bool: ``True`` se a sequência for não vazia e contiver apenas A/C/G/U;
        caso contrário ``False``.

    Raises:
        AttributeError: Se ``seq`` não for string e não tiver o método ``upper()``.

    Examples:
        >>> validar_rna("ACGU")
        True
        >>> validar_rna("ACGT")
        False
        >>> validar_rna("")
        False
    """
    if not seq:
        return False
    seq = seq.upper()
    return set(seq).issubset({'A', 'C', 'G', 'U'})


def validar_proteina(seq):
    """Valida se uma sequência contém apenas aminoácidos padrão (20 AA).

    Aceita o alfabeto:
    ``A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y`` (case-insensitive).
    Se a sequência for vazia (``""``/``None``/falsy), devolve ``False``.

    Args:
        seq (str): Sequência proteica a validar.

    Returns:
        bool: ``True`` se a sequência for não vazia e contiver apenas aminoácidos
        do conjunto aceite; caso contrário ``False``.

    Raises:
        AttributeError: Se ``seq`` não for string e não tiver o método ``upper()``.

    Examples:
        >>> validar_proteina("ACDEFGHIKLMNPQRSTVWY")
        True
        >>> validar_proteina("ACDXYZ")
        False
        >>> validar_proteina("")
        False
    """
    if not seq:
        return False
    seq = seq.upper()
    aminoacidos = {
        'A','C','D','E','F','G','H','I','K','L',
        'M','N','P','Q','R','S','T','V','W','Y'
    }
    return set(seq).issubset(aminoacidos)


def transcricao(dna_seq):
    """Transcreve DNA para RNA substituindo ``T`` por ``U``.

    Se a sequência de entrada não for DNA válido (ver :func:`validar_dna`),
    devolve string vazia ``""`` em vez de levantar exceção.

    Args:
        dna_seq (str): Sequência de DNA.

    Returns:
        str: Sequência de RNA resultante (em maiúsculas), ou ``""`` se a entrada
        não for DNA válido.

    Raises:
        AttributeError: Se ``dna_seq`` não for string e não tiver ``upper()``.

    Examples:
        >>> transcricao("ATGC")
        'AUGC'
        >>> transcricao("AUGC")
        ''
    """
    if not validar_dna(dna_seq):
        return ""
    return dna_seq.upper().replace('T', 'U')


def complemento(dna_seq):
    """Calcula o complemento de uma sequência de DNA.

    Usa o emparelhamento complementar:
    - A ↔ T
    - C ↔ G

    Se a sequência de entrada não for DNA válido (ver :func:`validar_dna`),
    devolve ``""``.

    Args:
        dna_seq (str): Sequência de DNA.

    Returns:
        str: Complemento (em maiúsculas), ou ``""`` se a entrada não for DNA válido.

    Raises:
        AttributeError: Se ``dna_seq`` não for string e não tiver ``upper()``.

    Examples:
        >>> complemento("ATGC")
        'TACG'
        >>> complemento("AUGC")
        ''
    """
    if not validar_dna(dna_seq):
        return ""
        
    mapa = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complemento = ""
    for base in dna_seq.upper():
        complemento += mapa[base]
    return complemento


def reverso(seq):
    """Devolve a sequência invertida (reverse).

    Se a sequência de entrada for vazia (``""``/``None``/falsy), devolve ``""``.

    Args:
        seq (str): Sequência a inverter.

    Returns:
        str: Sequência invertida, ou ``""`` se a entrada for vazia.

    Examples:
        >>> reverso("ATGC")
        'CGTA'
        >>> reverso("")
        ''
    """
    if not seq:
        return ""
    return seq[::-1]


def complemento_inverso(dna_seq):
    """Calcula o complemento inverso (reverse-complement) de uma sequência de DNA.

    Produz o complemento e inverte a orientação, equivalente ao complemento reverso
    usado frequentemente em bioinformática para obter a cadeia complementar no sentido
    5'→3'. Se a sequência de entrada não for DNA válido, devolve ``""``.

    Args:
        dna_seq (str): Sequência de DNA.

    Returns:
        str: Complemento inverso (em maiúsculas), ou ``""`` se a entrada não for DNA válido.

    Raises:
        AttributeError: Se ``dna_seq`` não for string e não tiver ``upper()``.

    Examples:
        >>> complemento_inverso("ATGC")
        'GCAT'
        >>> complemento_inverso("AUGC")
        ''
    """
    if not validar_dna(dna_seq):
        return ""
        
    mapa = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complemento = ""
    for base in dna_seq.upper():
        complemento += mapa[base]
    return complemento[::-1]


codao_stop = {"TAA", "TAG", "TGA"}

def encontra_codao_stop:
    for j in range(start_index, len(seq) -2,3):
        if seq[j:j+3] in codao_stop:
            return seq[start_index:j+3]

def get_orfs(dna):
    if not isinstance(dna, str):
        return []
    seq = dna.upper().replace(" ", "").replace("\n", "")
    orfs = []
    for i in range(0, len(seq) - 2,3):
        if seq[i:i+3] == "ATG":
            orf = encontra_codao_stop(seq,i)
            if orf:
                orfs.append(orf)
    return orfs
