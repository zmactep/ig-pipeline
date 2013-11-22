package common

import java.io.FileOutputStream

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 22.11.13
 * Time: 16:15
 */
object FileUtils {
  def writeFasta(out : FileOutputStream, data : Iterable[(String, String)]) : Unit = {
    data.foreach(rec => {
      val (seq, name) = rec
      out.write((">" + name + "\n").getBytes)
      out.write((seq + "\n\n").getBytes)
    })
  }
}
