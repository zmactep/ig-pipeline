import annotators.RegionAnnotator
import common.SequenceType
import igcont.ContainerUtils

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 9:59
 */
object Main {
  def main(args : Array[String]) = {
    val r = new RegionAnnotator("VDJH", if (args(1) == "amino") SequenceType.AMINO else SequenceType.NUCLEO, args(0))
    val query = args(2)
    var res : String = r.annotate(query).foldLeft("")((s, i) => s + i)
    println(res)
  }
}
