import alicont.algorithms.simple.LocalAlignment
import alicont.common.{Matrix, Scoring}
import common.Algo
import edu.msu.cme.rdp.readseq.readers.core.SFFCore
import edu.msu.cme.rdp.readseq.readers.{SequenceReader, Sequence}
import edu.msu.cme.rdp.readseq.writers.FastaWriter
import java.io.File
import scala.io.Source
import scala.util.{Failure, Success, Try}

class Splitter(midFile: File, scoringFile: File, debug: Boolean = false) {
  private val GAP_PENALTY = -3
  private val MIN_RESULT_LENGTH = 300

  private val mids = readMidFile(midFile)

  private val scoring = Scoring.loadMatrix(scoringFile.getName)

  def process_reads(input: File, outdir: File) = {
    def makeWriter(file: File) = {
      Try(new FastaWriter(file)) match {
        case Success(writer) =>
          writer
        case Failure(err) =>
          println(err.getMessage)
          sys.exit()
      }
    }

    var numSequences = 0

    // the following is a hack to get number of sequences in SFF file using internal API of ReadSeq library
    Try(new SFFCore(input)) match {
      case Success(core) =>
        numSequences = core.getCommonHeader.getNumReads
        core.close()
      case Failure(err) =>
        println(err.getMessage)
        sys.exit()
    }

    val reader = new SequenceReader(input) // now that we created the SFFCore, we may not check for the same exceptions again
    val processing_list = List.fill(numSequences)( () => reader.synchronized(reader.readNextSequence()) )
    val writers = mids.keys.map({key => (key, makeWriter(new File(outdir, s"$key.fasta")))}).toMap
    processing_list.par.foreach( extractor => process_read(extractor(), writers, debug))
    writers.foreach{case (_, w) => w.close()}
  }

  private def readMidFile(file: File) =
    Source.fromFile(file).getLines().map(_.split("\\s+").toList match {case id :: f :: r :: _ => (id, (new DnaString(f.trim), new DnaString(r.trim)))}).toMap

  private def getMidBounds(mid: String, midTrace: String, seqTrace: String) = {
    val midStart = midTrace.indexOf(mid.charAt(0))
    val midEnd = midTrace.lastIndexOf(mid.charAt(mid.length - 1)) + 1
    val gapsToStart = seqTrace.substring(0, midStart).count(_ == '-')
    val gapsToEnd = seqTrace.substring(0, midEnd).count(_ == '-')
    (midStart - gapsToStart, midEnd - gapsToEnd)
  }

  private def align(s: String, mid: String, m: Matrix) = {
    LocalAlignment.extendMatrix(mid, s, GAP_PENALTY, scoring, m)
    val (score, (seqTrace, midTrace)) = LocalAlignment.traceback(mid, s, GAP_PENALTY, scoring, m)
    val (midStart, midEnd) = getMidBounds(mid, midTrace, seqTrace)
    m.move(-mid.length)
    (score, midStart, midEnd, s"$midTrace\n$seqTrace")
  }

  private def alignMids(s: String, mid1: String, mid2: String, m: Matrix) =  {
    val (fScore, _, fEnd, debug1) = align(s, mid1, m)
    val (rScore, rStart, _, debug2) = align(s, mid2, m)

    (if (fEnd < rStart) fScore + rScore else -1.0, fEnd, rStart, s"$debug1\n$debug2")
  }

  private def process_read(sequence: Sequence, outputs: Map[String, FastaWriter], debug: Boolean = false) = {
    val seq = sequence.getSeqString
    val mat = new Matrix(seq.length, mids.map { case (_, (f, r)) => (f.length :: r.length :: Nil).max }.max + 1)
    val (((score,resultStart, resultEnd, debug_alignment), id), index) = mids.map({ case (id, (f, r)) =>
      (alignMids(seq, +f, -r, mat), id) :: (alignMids(seq, +r, -f, mat), id) :: Nil })
      .flatten
      .zipWithIndex
      .maxBy{ case (((s, _, _, _), _), _) => s }

    val valid = resultEnd - resultStart > MIN_RESULT_LENGTH
    val reversed = index % 2 == 1

    if (valid) {
      val result = new DnaString(seq.substring(resultStart, resultEnd))
      val output = outputs(id)
      output.synchronized {
        output.writeSeq(new Sequence(sequence.getSeqName, sequence.getDesc, if (reversed) -result else +result))
      }
    }

    if (debug) {
      synchronized {
        log(sequence.getSeqName, debug_alignment, score, resultStart, resultEnd, reversed, !valid)
      }
    }
  }

  private def log(name: String, debug_alignment: String, score: Double, resultStart: Int, resultEnd: Int, reversed:Boolean, skipped: Boolean) {
    val direction = if (reversed) "reverse" else "forward"
    val line = s">$name\n$debug_alignment\n$score $resultStart $resultEnd $direction"
    println(if (skipped) line.replaceAll("\n", "\nX\t") else line)
  }

  private class DnaString(seq: String) {
    def unary_+() = seq
    def unary_-() = Algo.reverseComp(seq)
    def length = seq.length
  }
}