package annotators

import igcont.{ContainerUtils, Container}
import common.{FileUtils, SequenceTrait}
import common.SequenceType.SequenceType
import alicont.AlignmentResult
import scala.collection.mutable.ArrayBuffer
import alicont.algorithms.AlgorithmType.AlgorithmType
import alicont.algorithms.AlgorithmType

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 12:02
 */
class RegionAnnotator(n : String, t : SequenceType,
                      gap : (Double, Double, Double) = (-10, -1, -5),
                      matrix : scala.Array[scala.Array[scala.Double]] = null,
                      algo : AlgorithmType = AlgorithmType.SEMIGLOBAL) {
  private val _regs     = Array("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")
  private val _name     = n
  private val _type     = new SequenceTrait(t)
  private val _cont     = new Container(_type.alphabet, _type.special, Array("Region"), _type.k)
  private val _algo     = algo
  private val _matrix   = matrix
  private val _gap_open = gap._1
  private val _gap_ext  = gap._2
  private val _gap_smpl = gap._3

  def this(n : String, t : SequenceType, fasta : String, kabat : String,
           gap : (Double, Double, Double) = (-10, -1, -5),
           matrix : scala.Array[scala.Array[scala.Double]] = null, algo : AlgorithmType = AlgorithmType.SEMIGLOBAL) = {
    this(n, t, gap, matrix, algo)
    _cont.addAnnotations("Region", _regs)

    // Load fasta
    FileUtils.readFasta(fasta).foreach(tpl => {
      val (name, seq) = tpl
      _cont.push(seq, name)
    })

    // Load kabat
    FileUtils.readKabat(kabat).foreach(tpl => {
      val (name, arr) = tpl
      val record = _cont.record(name)

      arr.zipWithIndex.foreach(tpl => {
        val ((start, end), reg) = tpl
        if (start != -1 && end != -1) {
          (start to end).foreach(pos => {
            record.setAnnotation(pos, 0, reg)
          })
        }
      })
    })
  }

  def find(pattern : String) : Iterable[(String, Int)] = _cont.find(pattern)

  def alignment(query : String, n : Int = 10) : Iterable[AlignmentResult] = {
    if (!AlgorithmType.affine.contains(_algo)) {
      _cont.alignment(query, _gap_smpl, if (_matrix == null) _type.score else matrix, _algo, n)
    }
    else {
      _cont.affine_alignment(query, _gap_open, _gap_ext, if (_matrix == null) _type.score else matrix, _algo, n)
    }
  }

  def annotate(query : String, n : Int = 3) : Iterable[Int] = {
    val anno = if (!AlgorithmType.affine.contains(_algo))
                 _cont.annotate(query, _gap_smpl, if (_matrix == null) _type.score else matrix, _algo, n)
               else
                 _cont.affine_annotate(query, _gap_open, _gap_ext, if (_matrix == null) _type.score else matrix, _algo, n)
    val result = ArrayBuffer.fill[Int](query.size)(0)

    anno.zipWithIndex.foreach(tpl => {
      val (node, i) = tpl
      if (node._2.contains("Region")) {
        result(i) = _regs.indexOf(node._2("Region"))
      }
      else {
        result(i) = 7
      }
    })

    result
  }

  def name : String = _name

  def stats() : Unit = ContainerUtils.print_stats(_cont)
}
