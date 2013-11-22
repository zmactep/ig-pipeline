package sangerrun

import org.biojava.bio.program.abi.ABITrace

import org.biojava3.alignment.{SimpleGapPenalty, Alignments, SubstitutionMatrixHelper}
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType
import org.biojava3.core.sequence.DNASequence

import java.io.{FileOutputStream, FilenameFilter, File}
import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer

import common.{FileUtils, Algo}


/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 21.11.13
 * Time: 10:12
 */
class PairSangerProcessing(input : String, output : String, primers : (String, String), local : Boolean = false) {
  private val _input   = input
  private val _output  = output
  private val _primers = primers
  private val _local   = local

  private val VL_LEADER_DEFAULT = "ATGAAATACCTGCTGCCGACCGCTGCTGCTGGTCTGCTGCTCCTCGCTGCCCAGCCGGCGATGGCTAGC"
  private val VH_LEADER_DEFAULT = "ATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCC"

  def process() : Unit = saveAll(processAll(getProcessingList))

  def getSequence(filename : String, rc : Boolean) : String = getSequence(new File(filename), rc)

  def getSequence(file : File, rc : Boolean) : String = {
    try {
      val abi = new ABITrace(file).getSequence.seqString().toUpperCase

      val r = if (!rc) abi else Algo.reverseComp(abi)
      if (!r.contains('N')) r else null
    }
    catch {
      case e : ArrayIndexOutOfBoundsException => return null
    }
  }

  def getFileList : mutable.WrappedArray[String] = {
    val dir = new File(_input)
    dir.list(new FilenameFilter {
      def accept(dir: File, name: String): Boolean = name.endsWith(".ab1")
    }).sorted
  }

  def getProcessingList : mutable.Iterable[(String, File, File)] = {
    val cache = Array[File](null, null)
    val processing_list = ArrayBuffer.empty[(String, File, File)]
    var cnt = 0
    getFileList.foreach(filename => {
      if (cnt % 2 == 0) {
        cache(0) = new File(_input, filename)
        cnt += 1
      }
      else {
        cache(1) = new File(_input, filename)
        val base = getBaseName(cache(0), cache(1))
        if (base == null) {
          cache(0) = cache(1)
          cache(1) = null
        }
        else {
          if (cache(0).getName.contains(_primers._1)) {
            processing_list += ((base, cache(0), cache(1)))
          }
          else {
            processing_list += ((base, cache(1), cache(0)))
          }
          cnt += 1
        }
      }
    })
    processing_list
  }

  def getBaseName(file1 : File, file2 : File) : String = {
    val s1 = file1.getName
    val s2 = file2.getName

    val s1pos = (s1.indexOf(_primers._1), s1.indexOf(_primers._2))
    val s2pos = (s2.indexOf(_primers._1), s2.indexOf(_primers._2))

    if (s1pos._1 != -1 && s1.substring(0, s1pos._1) == s2.substring(0, s2pos._2)) {
      s1.substring(0, s1pos._1)
    }
    else if (s1pos._2 != -1 && s1.substring(0, s1pos._2) == s2.substring(0, s2pos._1)) {
      s1.substring(0, s1pos._2)
    }
    else {
      null
    }
  }

  def saveAll(processed_all : mutable.Iterable[(String, String, String)]) : Unit = {
    val outfile = new File(_output)
    processed_all.foreach(tpl => {
      val (basename, vl, vh) = tpl
      val vl_n = (vl, basename + "-VL")
      val vh_n = (vh, basename + "-VH")

      var file = new File(outfile, basename + "-nucleo.fa")
      file.createNewFile()

      FileUtils.writeFasta(new FileOutputStream(file), Array(vl_n, vh_n))

      val vl_p = Algo.translateString(vl)
      val vh_p = Algo.translateString(vh)

      val vl_s = vl_p.indexOf('*')
      val vh_s = vh_p.indexOf('*')

      val vl_a = (vl_p.substring(0, if (vl_s != -1) vl_s else vl_p.size), basename + "-VL")
      val vh_a = (vh_p.substring(0, if (vh_s != -1) vh_s else vh_p.size), basename + "-VH")

      file = new File(outfile, basename + "-amino.fa")
      file.createNewFile()

      FileUtils.writeFasta(new FileOutputStream(file), Array(vl_a, vh_a))
    })
  }

  def processAll(all : mutable.Iterable[(String, File, File)]) : mutable.Iterable[(String, String, String)] =
    all.map(processOne).filter(_ != null)

  def processOne(one : (String, File, File)) : (String, String, String) = {
    val (basename, forward, backward) = one
    val forward_seq = getSequence(forward, rc = false)
    val backward_seq = getSequence(backward, rc = true)
    var sequence : String = null
    if (forward_seq == null && backward_seq != null) {
      sequence = backward_seq
    }
    else if (backward_seq == null && forward_seq != null) {
      sequence = forward_seq
    }
    else if (forward_seq == null || backward_seq == null) {
      return null
    }
    else {
      sequence = assemblePair(forward_seq, backward_seq)
    }

    val vlh = getVLH(sequence)
    if (vlh == null) {
      null
    }
    else {
      (basename, vlh._1, vlh._2)
    }
  }

  def getVLH(seq : String) : (String, String) = {
    var vl_start = 0
    var vl_end = 0
    var vh_start = 0
    if (_local) {
      val matrix = SubstitutionMatrixHelper.getNuc4_4
      val vl_align = Alignments.getPairwiseAlignment(new DNASequence(VL_LEADER_DEFAULT), new DNASequence(seq),
        PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(5, 2), matrix)
      val vh_align = Alignments.getPairwiseAlignment(new DNASequence(VH_LEADER_DEFAULT), new DNASequence(seq),
        PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(5, 2), matrix)

      vl_start = vl_align.getIndexInTargetAt(1) + vl_align.getSize - 1
      vl_end = vh_align.getIndexInTargetAt(1) - 1
      vh_start = vl_end + vh_align.getSize
    }
    else {
      val vl_lead_start = seq.indexOf(VL_LEADER_DEFAULT)
      if (vl_lead_start != -1) {
        vl_start = vl_lead_start + VL_LEADER_DEFAULT.length
        vl_end = seq.indexOf(VH_LEADER_DEFAULT)
        vh_start = vl_end + VH_LEADER_DEFAULT.length
      }
      else {
        vl_start = -1
        vl_end = -1
        vh_start = -1
      }
    }

    if (vl_start < 0 || vl_end < 0 || vh_start < 0) {
      null
    }
    else {
      (seq.substring(vl_start, vl_end), seq.substring(vh_start))
    }
  }

  def assemblePair(seq1 : String, seq2 : String) : String = {
    val (pos1, pos2, len) = Algo.longestSubstr(seq1, seq2)
    seq1.substring(0, pos1 + len) + seq2.substring(pos2 + len)
  }
}
