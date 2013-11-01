package common

import scala.io.Source

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 13:21
 */
object FileUtils {
  def readFasta(filename: String) : Stream[(String,String)] = {
    def parse(lines: Iterator[String]): Stream[(String,String)] = {
      if (lines.isEmpty) {
        return Stream[(String,String)]()
      }

      val name = lines.next().drop(1)
      val (seq,rest) = lines.span(_(0)!='>')
      (name, seq.mkString) #:: parse(rest)
    }
    parse(Source.fromFile(filename).getLines())
  }

  def readKabat(filename : String) : Stream[(String, Array[(Int, Int)])] = {
    def parse(lines: Iterator[String]) : Stream[(String, Array[(Int, Int)])] = {
      if (lines.isEmpty) {
        return Stream[(String, Array[(Int, Int)])]()
      }

      val jline = lines.next().mkString
      val line = jline.split('\t').filter(_ != "")
      val name = line(0)
      val arr = Array.fill[(Int, Int)](7)((0, 0))
      (0 to 6).foreach(i => arr(i) = (line(1 + 2*i).toInt - 1, line(2 + 2*i).toInt - 1))

      (name, arr) #:: parse(lines)
    }

    parse(Source.fromFile(filename).getLines())
  }
}
