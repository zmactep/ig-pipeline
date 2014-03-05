import java.io.File
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

      //just an example of output
      val seq: Sequence = "QLQLQESGGGLV"
      val triq = ProteinTriequence(seq.addAlternativeAt(0, Set("E")).addAlternativeAt(1, Set("V")).addAlternativeAt(4, Set("V")))
      println("Done constructing. Now printing to /tmp/graph.dot")
      triq.toDotFile(new File("/tmp/graph.dot"))
      println("Saved to file")
    }
  }
}
