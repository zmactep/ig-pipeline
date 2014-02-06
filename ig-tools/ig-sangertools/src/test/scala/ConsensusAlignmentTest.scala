import org.scalatest.FlatSpec
import org.scalatest.matchers.ShouldMatchers
import sangerrun.ConsensusAlignment

class ConsensusAlignmentTest extends FlatSpec with ShouldMatchers {
  "ConsensusAlignment" should "work on empty sequences" in {
    intercept[IllegalArgumentException] {ConsensusAlignment.merge() }
    ConsensusAlignment.merge("") should be ("")
    ConsensusAlignment.merge("", "") should be ("")
  }

  it should "align sequence on itself" in {
    ConsensusAlignment.merge("ABC") should be ("ABC")
    ConsensusAlignment.merge("ABC", "ABC") should be ("ABC")
    ConsensusAlignment.merge("ABC", "ABC", "ABC") should be ("ABC")
  }

  it should "align sequence on its subsequence" in {
    ConsensusAlignment.merge("SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQ", "YGFRLGFLHSGTAKS") should be ("SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQ")
    ConsensusAlignment.merge("LGFLH", "YGFRLGFLHSGTAKS", "SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQ") should be ("SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQ")
  }

  it should "merge reads" in {
    ConsensusAlignment.merge(
      "SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCS",
      "AIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEV",
      "DRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENL") should be
    ("SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENL")
  }
}
