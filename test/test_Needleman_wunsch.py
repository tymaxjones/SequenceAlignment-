from Needleman_wunsch import needleman_wunsch


def test_needleman_wunsch() -> None:
    """Test the needleman_wunsch function."""

    # Test case 1: Identical sequences
    assert needleman_wunsch("ACGTAT", "ACGTAT", 1, -1, -1) == (
        6,
        "ACGTAT",
        "ACGTAT",
    )

    # Test case 2: Different sequences
    assert needleman_wunsch("ACGTAT", "AGTGCT", 1, -1, -1) == (
        1,
        "ACGT-AT",
        "A-GTGCT",
    )

    # Test case 3: Unequal length sequences
    assert needleman_wunsch("ACGTAT", "ACG", 1, -1, -1) == (
        0,
        "ACGTAT",
        "ACG---",
    )

    # Test case 4: Different gap penalties
    assert needleman_wunsch("ACGTAT", "AGTGCT", 1, -1, -1) == (
        1,
        "ACGT-AT",
        "A-GTGCT",
    )

    # Test case 5: Verbose Output
    assert needleman_wunsch("ACTGC", "ACTCA", 1, -1, -1, True) == (
        "Alignment Score: 2\nACTGC-\n||| | \nACT-CA"
    )


test_needleman_wunsch()
