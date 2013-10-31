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
    ContainerUtils.warmup()

    println("Start!")
    val start = System.currentTimeMillis()
    //val a = new Annotator("VDJH", SequenceType.NUCLEO, "../../data/train/VDJH_train")
    val a = new Annotator("VDJH", SequenceType.NUCLEO, "../../data/BIG_FILE")
    printf("Time: %f", (System.currentTimeMillis() - start) / 1000.0)

    a.stats()

    println(a.name)
  }
}
