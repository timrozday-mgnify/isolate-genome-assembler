from pathlib import Path


def test_empty_samplesheet_is_explicit() -> None:
    """The Phase 0 stub test must have a valid, intentionally empty input."""
    samplesheet = Path(__file__).parent / "data" / "empty_samplesheet.yml"
    assert samplesheet.read_text(encoding="utf-8") == "samples: []\n"


def test_hifi_fixture_has_a_complete_fastq_record() -> None:
    reads = Path(__file__).parent / "data" / "isolate01.hifi_reads.fastq"
    assert reads.read_text(encoding="utf-8").splitlines() == [
        "@read_001",
        "ACGT",
        "+",
        "####",
    ]
