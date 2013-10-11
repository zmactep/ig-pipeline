package utils

import java.io.File

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 11.10.13
 * Time: 9:12
 * To change this template use File | Settings | File Templates.
 */

object FileUtils {
  def fixPath(path: String) : String = if (path.endsWith("/")) path else path + "/"

  def createDirIfNotExists(path: String) = {
    val dirExistCheck = new File(path)
    if (!dirExistCheck.exists()) {
      dirExistCheck.mkdir()
    }
  }
}
