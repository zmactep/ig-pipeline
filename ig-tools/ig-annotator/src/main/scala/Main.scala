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
    ContainerUtils.warmup()

    val r = new RegionAnnotator("VDJH", SequenceType.NUCLEO, "../../data/train/VDJH_train")
    r.stats()

    val query = args(0)
    var time = System.currentTimeMillis()
    var res : String = null
    res = r.annotate(query).foldLeft("")((s, i) => s + i)
    println(query)
    println(res)
  }
}
