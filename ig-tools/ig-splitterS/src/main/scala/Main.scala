import java.io.{FileNotFoundException, File}
import scala.util.{Failure, Try}
import scopt.OptionParser

object Main {
  private val MATRIX_FILENAME = "NUC1.1.txt"
  private val MIDS_FILENAME = "mids.txt"

  private case class Config(input: File = null, outdir: File = null, mids: File = new File(MIDS_FILENAME), scoring: File = new File(MATRIX_FILENAME), debug: Boolean = false)

  private def getParser : scopt.OptionParser[Config] = new OptionParser[Config]("ig-splitterS") {
    head("ig-splitterS", "1.0-SNAPSHOT")
    note("Required:")
    arg[File]("input") action { (x, c) => c.copy(input = x) } text "input sff filename to split"
    arg[File]("outdir") action { (x, c) => c.copy(outdir = x) } validate { x => if (x.isDirectory) success else failure("Not a directory")} text "output directory"
    note("Options:")
    opt[File]("mids") action { (x, c) => c.copy(mids = x) } text s"mids filename (default=$MIDS_FILENAME)"
    opt[File]("scoring_matrix") action { (x, c) => c.copy(scoring = x) } text s"scoring matrix filename (default=$MATRIX_FILENAME)"
    opt[Unit]("debug") action { (_, c) => c.copy(debug = true) } text "enable debug output"
    help("help") text "this message"
  }

  def main(args: Array[String]) = {
    val parser = getParser

    parser.parse(args, Config()) map {config => {
      Try( new Splitter(config.mids, config.scoring, config.debug).process_reads(config.input, config.outdir) ) match {
        case Failure(e: FileNotFoundException) =>
          Console.err.printf("One of input files was not found: %s\n", e.getMessage)
          sys.exit(0)
        case Failure(e : Exception) =>
          Console.err.printf("Unknown fatal error: %s\n", e)
          sys.exit(0)
        case _ => ()
      }
    }} getOrElse {
      parser.showUsage
    }
  }
}

