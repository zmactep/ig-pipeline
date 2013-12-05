import annotators.RegionAnnotator
import common.SequenceType
import java.io.File

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 9:59
 */
object Main {
  case class Config(amino : Boolean = false, fasta : String = null, kabat : String = null, source : String = null)

  def main(args : Array[String]) = {
    val parser = new scopt.OptionParser[Config]("ig-annotator") {
      head("ig-annotator", "1.0-SNAPSHOT")
      opt[String]('s', "source") required() action {(s, c) => c.copy(source = s)} text("")
      opt[Unit]("amino") action {(_, c) => c.copy(amino = true)} text("use amino acid data")
    }

    val r = new RegionAnnotator("Name", if (args(1) == "amino") SequenceType.AMINO else SequenceType.NUCLEO, args(0))
    val query = args(2)
    val res : String = r.annotate(query).foldLeft("")((s, i) => s + i)
    println(res)
  }
}
