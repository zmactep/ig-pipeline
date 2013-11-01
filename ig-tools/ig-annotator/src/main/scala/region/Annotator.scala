package region

import igcont.{ContainerUtils, Container}
import common.{FileUtils, SequenceTrait}
import common.SequenceType.SequenceType
import alicont.AlignmentResult
import igcont.anno.Record
import scala.collection.immutable.HashMap

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 12:02
 */
class Annotator(n : String, t : SequenceType) {
  private val _regs = Array("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")
  private val _name = n
  private val _type = new SequenceTrait(t)
  private val _cont = new Container(_type.alphabet, _type.special, Array("Region"), _type.k)

  def this(n : String, t : SequenceType, fileprefix : String) = {
    this(n, t)
    _cont.addAnnotations("Region", _regs)

    // Load fasta
    FileUtils.readFasta(fileprefix + ".fasta").foreach(tpl => {
      val (name, seq) = tpl
      _cont.push(seq, name)
    })

    // Load kabat
    FileUtils.readKabat(fileprefix + ".kabat").foreach(tpl => {
      val (name, arr) = tpl
      val record = _cont.record(name)

      arr.zipWithIndex.foreach(tpl => {
        val ((start, end), reg) = tpl
        (start to end).foreach(pos => {
          record.setAnnotation(pos, 0, reg)
        })
      })
    })
  }

  def find(pattern : String) : Iterable[(String, Int)] = _cont.find(pattern)

  def alignment(query : String, n : Int = 10) : Iterable[AlignmentResult] = _cont.alignment(query, -5, _type.score, n)

  def name : String = _name

  def stats() : Unit = ContainerUtils.print_stats(_cont)
}