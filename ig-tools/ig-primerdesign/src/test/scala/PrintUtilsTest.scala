import org.scalatest.FunSpec
import ru.biocad.ig.primer._
import ru.biocad.ig.primer.ProteinTriequence

/**
 * @author kfeodorov
 * @since 22.02.14
 */
class PrintUtilsTest extends FunSpec {
  describe("PrintUtils") {
    it ("should work with Pasha's dataset :)") {
      import ru.biocad.ig.primer.ConversionUtils._
      import ru.biocad.ig.primer.DnaUtils._
      for (protein <- List("DIQMTQSPSSLSASVGDRVTITCKASQSVSSDVGWYQQKPGKAPKLLIYSGSNRYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQDYYSPWTFGQGTKVEIK","EVQLVQSGAEVKKPGSSVKVSCKASGYTFTNYVINWVRQAPGQGLEWIGYNDPYNDVSKYNEKFKGRATITSDKSTSTAYMELSSLRSEDTAVYYCAKEGGGKYVYAMDSWGQGTTVTVSS")) {
        val triq = ProteinTriequence(protein)
        val ds = new SimpleProbabilityDecisionStrategy(Map("A" -> 5, "C" -> 5, "G" -> 5, "T" -> 5, "U" -> 5))
        println(protein)
        val nucl: Option[Sequence] = triq.sample(ds).map{s => println(s"$s straight"); println(s"${reverseComplementRNA(Option(s)).get} revcomp"); String2Sequence(s)}
        val indicies = DnaUtils.findOverlaps(nucl, 20, 40)
        PrintUtils.printOverlaps(nucl, indicies, 20)
      }
    }
  }
}