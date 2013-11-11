package utils

import java.io._
import com.typesafe.config.ConfigFactory
import scala.sys.process.Process
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
      dirExistCheck.mkdirs
    }
  }

  def getBiggestNameDir(path: String): String = {
    val dir = new File(path)
    if (! dir.exists()) {
      throw new FileNotFoundException
    }
    
    val dirs = dir.listFiles(
      new FilenameFilter() {
        @Override
        override def accept(dir: File, name: java.lang.String): Boolean = {
          name.matches("\\d+")
        }
      }
    ).toList.filter(_.isDirectory)
    if (dirs.length == 0) {
      throw new FileNotFoundException
    }

    var lastModifiedDir = dirs(0)
    dirs.drop(0).foreach { f =>
      if (lastModifiedDir.getName.toInt < f.getName.toInt) {
        lastModifiedDir = f
      }
    }

    lastModifiedDir.getName
  }

  def saveToFile(path: String, data: String) = {
    try {
      val out = new BufferedWriter(new FileWriter(new File(path)))
      out.write(data)
      out.close()
    } catch {
      case e: IOException =>
    }
  }

  def scanDir(path: String, toolsRoot: String, group: String, run: String) = {
    val conf = ConfigFactory.load()
    val params = s"--host=${conf.getString("ig-backend.db_host")} --port=${conf.getString("ig-backend.db_port")} " +
      s"--db=${conf.getString("ig-backend.db_name")} --user=${conf.getString("ig-backend.db_user")} " +
      s"--password=${conf.getString("ig-backend.db_password")} --dir=$path --group=$group --run=$run"

    val cmd = s"python ${new File(toolsRoot, "ig-snooper/ig_snooper_utils/scan_dir.py").toString} $params"
    Process(cmd, new java.io.File(toolsRoot)).!!
  }
}
