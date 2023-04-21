from Smith_waterman import smith_waterman


def test_smith_waterman() -> None:
    """Test the smith_waterman function."""

    # Test case 1: Identical sequences
    assert smith_waterman("ACGTAT", "ACGTAT", 10, -5, -5) == (
        100,
        60,
        "ACGTAT",
        "ACGTAT",
        "ACGTAT",
    )

    # Test case 2: Different sequences
    assert smith_waterman("ACGTAT", "AGTGCT", 10, -5, -5) == (
        57,
        25,
        "ACGT-AT",
        "A GT  T",
        "A-GTGCT",
    )

    # Test case 3: Unequal length sequences
    assert smith_waterman("ACGTAT", "ACG", 10, -5, -5) == (
        100,
        30,
        "ACG",
        "ACG",
        "ACG",
    )

    # Test case 4: Different gap penalties
    assert smith_waterman("ACGTAT", "AGTGCT", 10, -5, -3) == (
        57,
        29,
        "ACGT-AT",
        "A GT  T",
        "A-GTGCT",
    )


test_smith_waterman()
