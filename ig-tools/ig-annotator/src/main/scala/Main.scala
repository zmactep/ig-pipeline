import alicont.algorithms.AlgorithmType._
import alicont.common.Scoring
import annotators.RegionAnnotator
import common.{FileUtils, SequenceType}
import java.io.FileNotFoundException

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 9:59
 */
object Main {
  case class Config(amino : Boolean = false, fasta : String = null, kabat : String = null, source : String = null,
                    gap_open : Double = -10, gap_ext : Double = -1, gap : Double = -5,
                    matrix : String = null, align : AlgorithmType = SEMIGLOBAL, marking : Boolean = false,
                    filter : Int = 0)

  def getParser : scopt.OptionParser[Config] = new scopt.OptionParser[Config]("ig-annotator") {
    head("ig-annotator", "1.0-SNAPSHOT")
    note("Required:")
    opt[String]('s', "source") required() action {(s, c) => c.copy(source = s)} text "file to annotate [fasta]"
    opt[String]('r', "reference") required() action {(s, c) => c.copy(fasta = s)} text "reference file [fasta]"
    opt[String]('m', "marking") required() action {(s, c) => c.copy(kabat = s)} text "reference marking [igblast marking format]"
    note("Optional:")
    opt[Unit]("amino") action {(_, c) => c.copy(amino = true)} text "use amino acid data"
    opt[Unit]("igblast-like") action {(_, c) => c.copy(marking = true)} text "output as igblast marking"
    opt[Int]("filter-window") action {(s, c) => c.copy(filter = s)} text "set window size for filtration (default: disabled)"
    note("  alignment parameters")
    opt[String]("matrix") action {(s, c) => c.copy(matrix = s)} text "use external alignment matrix [txt]"
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

  def constructAnnotator(config : Config) : RegionAnnotator =
    new RegionAnnotator("Name", if (config.amino) SequenceType.AMINO else SequenceType.NUCLEO,
                        config.fasta, config.kabat,
                        (config.gap_open, config.gap_ext, config.gap),
                        if (config.matrix != null) Scoring.loadMatrix(config.matrix) else null,
                        config.align)

  def constructMarkup(anno : String) : String = {
    val r = Array.fill[Array[Int]](7)(Array.fill[Int](2)(0))
    val sb = new StringBuilder()
    anno.zipWithIndex.foreach(tpl => {
      val (c, i) = tpl
      val ci = c.asDigit
      if (ci != 7) {
        if (r(ci)(0) == 0) {
          r(ci)(0) = i + 1
        }
        else {
          r(ci)(1) = i + 1
        }
      }
    })
    (0 until r.size).foreach(i => {
      sb.append("\t%d\t%d".format(r(i)(0), r(i)(1)))
    })

    sb.toString()
  }

  def run(config : Config) : Unit = {
    val r = constructAnnotator(config)
    val data = FileUtils.readFasta(config.source)
    data.foreach(tpl => {
      val (name, seq) = tpl
      val raw_anno = r.annotate(seq).foldLeft("")((s, i) => s + i)
      val anno = if (config.filter > 1) filterAnno(raw_anno, config.filter) else raw_anno

      if (!config.marking) {
        printf("> %s\n%s\n%s\n\n", name, seq, anno)
      }
      else {
        printf("%s%s\n", name, constructMarkup(anno))
      }
    })
  }

  def filterAnno(anno : String, window : Int) : String = {
    val sb = new StringBuilder()
    (0 until anno.size).foreach(i => {
      sb += anno.substring(math.max(0, i - window / 2),
                           math.min(anno.size, i + window / 2)).groupBy(identity).maxBy(_._2.size)._1
    })
    sb.toString()
  }

  def main(args : Array[String]) = {
    val parser = getParser

    parser.parse(args, Config()) map {config => {
      try {
        run(config)
      } catch {
        case e : FileNotFoundException => {
          Console.err.printf("One of input files was not found: %s\n", e.getMessage)
          sys.exit(0)
        }
        case e : Exception => {
          Console.err.printf("Unknown fatal error: %s\n", e)
          sys.exit(0)
        }
      }
    }} getOrElse {
      parser.showUsage
    }
  }
}
