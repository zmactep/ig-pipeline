import java.io.File
import java.util.Collections
import org.apache.commons.io.FileUtils
import org.scalatest.FunSpec
import org.scalatest.Matchers._
import ru.biocad.ig.primer.DnaUtils._
import ru.biocad.ig.primer._
import ru.biocad.ig.primer.ProteinTriequence

/**
 * @author kfeodorov
 * @since 02.03.14
 */
class TriquenceTest extends FunSpec {
  describe("Triquence") {
    import ru.biocad.ig.primer.ConversionUtils._

    it ("correctly handle input") {
      Triequence(Option("EW")).asInstanceOf[ProteinTriequence].protein should be (List(Set("E"), Set("W")))
      Triequence(None) should be (EmptyTriquence)

      val seq: Sequence = "EW"
      val triq: Triequence = Triequence(Option(seq.addAlternativeAt(0, Set("S"))))
      triq.asInstanceOf[ProteinTriequence].protein should be (List(Set("E", "S"), Set("W")))
    }

    it ("should build an empty Triequence") {
      ProteinTriequence("EW").protein should be (List(Set("E"), Set("W")))
      Triequence(None) should be (EmptyTriquence)
    }

    val protein: String = "QLQLQESGGGLVQPGGSLRLSCAASGFTFTRYAMSWVRQAPGKGLEWVSAINSGGGSTYYADSVKGRFTISRDNAKNTVYLQLNSLKTEDMADVLVCRGGGTLGGR"
    val seq: Sequence = protein

    it ("should build a Triequence") {
      //just an example of output. Sequence is taken from readme
      val triq = ProteinTriequence(seq.addAlternativeAt(0, Set("E")).addAlternativeAt(1, Set("V")).addAlternativeAt(4, Set("V")).
        addAlternativeAt(29, Set("S")).addAlternativeAt(74, Set("S")).addAlternativeAt(78, Set("L")).addAlternativeAt(82, Set("M")).
        addAlternativeAt(86, Set("R")).addAlternativeAt(87, Set("A")).addAlternativeAt(90, Set("T")).addAlternativeAt(94, Set("Y")).
        addAlternativeAt(95, Set("Y")))
      FileUtils.writeLines(new File("/tmp/graph.dot"), Collections.singletonList(triq), false)
    }

    it ("should correctly sample nucleotide strings") {
      val ds = new SimpleProbabilityDecisionStrategy(Map("A" -> 5, "C" -> 5, "G" -> 5, "T" -> 5, "U" -> 5))
      val triq = ProteinTriequence(seq)
      def translate(rna: String) = augmentString(rna).sliding(3, 3).toList.map{codon2aminoAcid}.mkString("")
      (1 to 1000).foreach(_ => translate(triq.sample(ds).get) should be (protein))
    }
  }
}
