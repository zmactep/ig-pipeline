import sangerrun.PairSangerProcessing

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 21.11.13
 * Time: 9:55
 */
object Main {
  def main(args : Array[String]) = {
    val p = new PairSangerProcessing("/home/mactep/Data/abi-test", "/home/mactep/Data/abi-test/result", ("pLLSeqR", "H3b"), true)
    p.process()
  }
}
