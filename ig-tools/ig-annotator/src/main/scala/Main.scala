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
  def main(s : Array[String]) = {
    ContainerUtils.warmup()

    val r = new RegionAnnotator("VDJH", SequenceType.NUCLEO, "../../data/train/VDJH_train")
    r.stats()

    val query = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCGGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGACGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCAGTACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGTCGTGTATTACTGTGCGAGAGATAGTAGCTCCCACTATACCCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG"
    var time = System.currentTimeMillis()
    var res : String = null
    (1 to 100).foreach(i => {
      res = r.annotate(query).foldLeft("")((s, i) => s + i)
    })
    time = System.currentTimeMillis() - time
    printf("Time: %.2f\n", time / 1000.0)
    printf("Avg time: %.2f\n", time / 1000.0 / 100)
  }
}
