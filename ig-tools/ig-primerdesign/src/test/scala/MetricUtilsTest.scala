import org.scalatest.FunSpec
import org.scalatest.Matchers._
import ru.biocad.ig.primer.MetricUtils

/**
 * @author kfeodorov
 * @since 21.03.14
 */
class MetricUtilsTest extends FunSpec {
  describe("MetricUtils") {
    it("correctly calculates score") {
      MetricUtils.similarityScore(Option("AATT"), Option("AATT")) should be (Some((8, "AATT")))
      MetricUtils.similarityScore(Option("CCGG"), Option("CCGG")) should be (Some((12, "CCGG")))
      MetricUtils.similarityScore(Option("AAAA"), Option("TTTT")) should be (Some((0,"")))
      MetricUtils.similarityScore(Option("ACCGGTTT"), Option("CGGT")) should be (Some((11,"CGGT")))
    }

    it("correctly calculates hairpin score as in Pasha's toolkit") {
      MetricUtils.hairpinScore(Option("AAAATTTT")) should be (Some(1.0))
      MetricUtils.hairpinScore(Option("AAAACTTTT")) should be (Some (0.42))
      MetricUtils.hairpinScore(Option("AAGATTCTT")) should be (9)
      MetricUtils.hairpinScore(Option("AAGGGTCTT")) should be (7)
    }
  }
}