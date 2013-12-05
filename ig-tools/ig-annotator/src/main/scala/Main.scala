import alicont.algorithms.AlgorithmType._
import alicont.common.Scoring
import annotators.RegionAnnotator
import common.{FileUtils, SequenceType}
import java.io.{FileNotFoundException, File}

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 9:59
 */
object Main {
  case class Config(amino : Boolean = false, fasta : String = null, kabat : String = null, source : String = null,
                    gap_open : Double = -10, gap_ext : Double = -1, gap : Double = -5,
                    matrix : String = null, align : AlgorithmType = SEMIGLOBAL)

  def main(args : Array[String]) = {
    val parser = new scopt.OptionParser[Config]("ig-annotator") {
      head("ig-annotator", "1.0-SNAPSHOT")
      note("Required:")
      opt[String]('s', "source") required() action {(s, c) => c.copy(source = s)} text "fasta file to annotate"
      opt[String]('r', "reference") required() action {(s, c) => c.copy(fasta = s)} text "reference fasta file"
      opt[String]('m', "marking") required() action {(s, c) => c.copy(kabat = s)} text "reference marking"
      note("Optional:")
      opt[Unit]("amino") action {(_, c) => c.copy(amino = true)} text "use amino acid data"
      note("  alignment parameters")
      opt[String]("matrix") action {(s, c) => c.copy(matrix = s)} text "use external alignment matrix"
      opt[Double]("gap") action {(s, c) => c.copy(gap = s)} text "simple gap score (default: -5)"
      opt[Double]("gap-open") action {(s, c) => c.copy(gap_open = s)} text "affine open gap score (default: -10)"
      opt[Double]("gap-ext") action {(s, c) => c.copy(gap_ext = s)} text "affine extension gap score (default: -1)"
      note("  alignment algorithms")
      opt[Unit]("global") action {(_, c) => c.copy(align = GLOBAL)} text "use global alignment"
      opt[Unit]("local") action {(_, c) => c.copy(align = LOCAL)} text "use local alignment"
      opt[Unit]("semiglobal") action {(_, c) => c.copy(align = SEMIGLOBAL)} text "use semiglobal alignment (default)"
      opt[Unit]("affine-global") action {(_, c) => c.copy(align = AFFINE_GLOBAL)} text "use global alignment"
      opt[Unit]("affine-local") action {(_, c) => c.copy(align = AFFINE_LOCAL)} text "use global alignment"
      opt[Unit]("affine-semiglobal") action {(_, c) => c.copy(align = AFFINE_SEMIGLOBAL)} text "use global alignment"
      help("help") text "this message"
    }

    parser.parse(args, Config()) map {config => {
      try {
        val matrix = if (config.matrix != null) Scoring.loadMatrix(config.matrix) else null
        val r = new RegionAnnotator("Name", if (config.amino) SequenceType.AMINO else SequenceType.NUCLEO,
                                    config.fasta, config.kabat,
                                    (config.gap_open, config.gap_ext, config.gap),
                                    matrix, config.align)
        val data = FileUtils.readFasta(config.source)
        data.foreach(tpl => {
          val (name, seq) = tpl
          val anno = r.annotate(seq).foldLeft("")((s, i) => s + i)

          printf("> %s\n%s\n%s\n\n", name, seq, anno)
        })
      } catch {
        case e : FileNotFoundException => {
          Console.err.printf("One of input files was not found: %s\n", e.getMessage)
          sys.exit(0)
        }
        case e : Exception => {
          Console.err.printf("Unknown fatal error: %s\n", e.getMessage)
          sys.exit(0)
        }
      }
    }} getOrElse {
      parser.showUsage
    }
  }
}
