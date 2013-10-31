package region

import igcont.{ContainerUtils, Container}
import common.{FileUtils, SequenceTrait}
import common.SequenceType.SequenceType

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 12:02
 */
class Annotator(n : String, t : SequenceType) {
  private val _name = n
  private val _type = new SequenceTrait(t)
  private val _cont = new Container(_type.alphabet, _type.special, Array("Region"), _type.k)

  def this(n : String, t : SequenceType, fileprefix : String) = {
    this(n, t)
    _cont.addAnnotations("Region", Array("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"))

    val rc = Runtime.getRuntime

    FileUtils.readFasta(fileprefix + ".fasta").takeWhile(_ => (rc.totalMemory() - rc.freeMemory()) / 1024 / 1024 < 1500).foreach(tpl => {
      val (name, seq) = tpl
      _cont.push(seq, name)
    })
  }

  def name : String = _name

  def stats() : Unit = ContainerUtils.print_stats(_cont)
}
