package common

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 13:21
 */
object FileUtils {
  def readFasta(filename: String) = {
    import scala.io.Source
    def parse(lines: Iterator[String]): Stream[(String,String)] = {
      if(lines.isEmpty) return Stream[(String,String)]()
      val name = lines.next().drop(1)
      val (seq,rest) = lines.span(_(0)!='>')
      (name, seq.mkString) #:: parse(rest)
    }
    parse(Source.fromFile(filename).getLines())
  }
}
