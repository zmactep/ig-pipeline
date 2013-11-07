package workers

import sys.process._
import scala.concurrent.Future
import akka.actor.{ActorPath, ActorRef, ActorLogging}
import akka.pattern.pipe
import com.googlecode.protobuf.format.JsonFormat
import java.io._
import utils.FileUtils
import protocol.Command.{ResponseCommand, BatchCommand}
import scala.collection.JavaConverters._
import protocol.Command.BatchCommand.RequestCommand


/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 25.09.13
 * Time: 10:02
 * To change this template use File | Settings | File Templates.
 */

class SimpleWorker(masterLocation: ActorPath) extends Worker(masterLocation) with ActorLogging{
  implicit val ec = context.dispatcher
  val toolsRoot = FileUtils.fixPath(context.system.settings.config.getString("ig-backend.tools_root"))
  val storageRoot = FileUtils.fixPath(context.system.settings.config.getString("ig-backend.storage_root"))

  FileUtils.createDirIfNotExists(storageRoot)

  def doWork(workSender: ActorRef, msg: Any): Unit = {
    Future {
      val batchBuilder = BatchCommand.newBuilder()
      try {
        msg match {
          case (m, _) => {
            JsonFormat.merge(m.toString, batchBuilder)
            val results = batchBuilder.build.getCommandsList.asScala.map(executeCommand).mkString("\n")
            self ! WorkComplete(results)
          }
        }
      } catch {
        case e: Exception => {
          log.debug(s"Error for task ${msg}: ${e.toString}")
          self ! akka.actor.Status.Failure(e)
        }
      }
    } pipeTo self
  }

  private def executeCommand(command: RequestCommand): String = {
    def buildCommand(outDir: String): String = {
      val executable = new File(toolsRoot, command.getExecutable).toString
      val params = command.getInput.getParamsList.asScala.map(p => s"--${p.getName}=${p.getValue}").mkString(" ")

      val cmd = s"python $executable $params --tools_root=${toolsRoot} --outdir=$outDir"
      log.debug(cmd)
      cmd
    }

    def formatOutput(response: String) : String = {
      val responseBuilder = ResponseCommand.newBuilder()
      //TODO use better way to check if everything is OK
      if (response.contains("Done in ")) responseBuilder.setStatus("ok") else responseBuilder.setStatus("failed")
      responseBuilder.setMessage(response)
      JsonFormat.printToString(responseBuilder.build())
    }

    val group = command.getInput.getGroup
    //runNumber == subdir name in group. Name is just  an increasing number. So, latestModified dir == dir with largest name
    var runNumber = 0

    try {
      runNumber = FileUtils.getLastModifiedDir(new File(storageRoot, group).toString).toInt + 1
    } catch {
      case _: NumberFormatException => log.debug("No numeric folders in " + new File(storageRoot, group).toString)
      case _: FileNotFoundException => log.debug("No last modified folder in " + new File(storageRoot, group).toString)
    }

    val outDir = new File(storageRoot, new File(group, runNumber.toString).toString).toString
    FileUtils.createDirIfNotExists(outDir)
    val execResult = Process(buildCommand(outDir), new java.io.File(toolsRoot)).!!
    FileUtils.scanDir(outDir, toolsRoot, group, runNumber.toString)

    FileUtils.saveToFile(new File(outDir, "description.txt").toString, execResult)
    val response = formatOutput(execResult)
    log.debug("Response: " + response)
    response
  }
}

