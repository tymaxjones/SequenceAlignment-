from Multiple_alignment import *


def test_Multiple_alignment() -> None:
    """Test the Multiple_alignment function."""

    # Test case 1: Identical sequences
    assert multiple_alignment(["ATCG", "CTCG", "GCAT", "ATCG"]) == [
        "ATCG",
        "ATCG",
        "CTCG",
        "AT--",
    ]

    # Test case 2: Different sequences
    assert multiple_alignment(["GCAT", "ATCG"]) == ["GCAT--", "--ATCG"]

    # Test case 3: Unequal length sequences
    assert multiple_alignment(["GCAT", "ATCG", "CATG"]) == [
        "GCAT-",
        "-CATG",
        "-ATCG",
    ]

    # Test case 5: Verbose Output
    assert multiple_alignment(["GCAT", "ATCG", "CATG", "ATCG", "ATCG"]) == [
        "ATCG",
        "ATCG",
        "ATCG",
        "AT-G",
        "AT--",
    ]


test_Multiple_alignment()
