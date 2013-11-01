import igcont.ContainerUtils
import region.Annotator
import common.SequenceType

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 9:59
 */
object Main {
  def main(s : Array[String]) = {
    println("Warmup")
    //ContainerUtils.warmup()

    println("Start!")
    val start = System.currentTimeMillis()
    val a = new Annotator("VDJH", SequenceType.NUCLEO, "../../data/train/VDJH_train")

    printf("Time: %f\n", (System.currentTimeMillis() - start) / 1000.0)
    a.stats()

    val res = a.alignment("CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCGGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGACGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCAGTACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGTCGTGTATTACTGTGCGAGAGATAGTAGCTCCCACTATACCCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG")

    res.foreach(aln => printf("\nName: %s\nScore: %d (%.2f)\nQ: %s\nT: %s\n", aln.name, aln.score, aln.similarity, aln.query, aln.target))

    println(a.name)
  }
}
