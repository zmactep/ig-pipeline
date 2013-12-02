import sangerrun.PairSangerProcessing

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 21.11.13
 * Time: 9:55
 */
object Main {
  def main(args : Array[String]) = {
    val p = new PairSangerProcessing("/Users/pavel/BIO/TARGETS/CMET/CMETHER3_MMMP1-raw data all/RAW",
                                     "/Users/pavel/BIO/TARGETS/CMET/CMETHER3_MMMP1-raw data all/RAW/result",
                                     ("SeqR", "H3b"), true)
    p.process()
  }
}
