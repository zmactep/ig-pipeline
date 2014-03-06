import java.io.File
import java.util.Collections
import org.apache.commons.io.FileUtils
import org.scalatest.FunSpec
import org.scalatest.Matchers._
import ru.biocad.ig.primer.DnaUtils._
import ru.biocad.ig.primer.{Triequence, EmptyTriquence, ProteinTriequence}

/**
 * Created by Kos on 02.03.14.
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

    it ("should build a trie") {
      ProteinTriequence("EW").protein should be (List(Set("E"), Set("W")))
      Triequence(None) should be (EmptyTriquence)

      //just an example of output. Sequence is taken from readme
      val seq: Sequence = "QLQLQESGGGLVQPGGSLRLSCAASGFTFTRYAMSWVRQAPGKGLEWVSAINSGGGSTYYADSVKGRFTISRDNAKNTVYLQLNSLKTEDMADVLVCRGGGTLGGR"
      val triq = ProteinTriequence(seq.addAlternativeAt(0, Set("E")).addAlternativeAt(1, Set("V")).addAlternativeAt(4, Set("V")).
        addAlternativeAt(29, Set("S")).addAlternativeAt(74, Set("S")).addAlternativeAt(78, Set("L")).addAlternativeAt(81, Set("M")).
        addAlternativeAt(84, Set("R")).addAlternativeAt(85, Set("A")).addAlternativeAt(88, Set("T")).addAlternativeAt(92, Set("Y")).
        addAlternativeAt(93, Set("Y")))
      FileUtils.writeLines(new File("/tmp/graph.dot"), Collections.singletonList(triq), false)
    }
  }
}
