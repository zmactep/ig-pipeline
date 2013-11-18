package workers

import sys.process._
import scala.concurrent.Future
import akka.actor.{ActorPath, ActorRef, ActorLogging}
import akka.pattern.pipe
import com.googlecode.protobuf.format.JsonFormat
import java.io._
import utils.FileUtils
import protocol.Command.{ResponseBatch, BatchCommand}
import protocol.Command.BatchCommand.RequestCommand
import scala.collection.JavaConverters._
import scala.collection.mutable.ListBuffer
import protocol.Command.ResponseBatch.ResponseCommand
import protocol.Command.ResponseBatch.ResponseCommand.{Builder, Files}


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
            val responseBatch = ResponseBatch.newBuilder()
            val outputDirs = ListBuffer[String]()
            batchBuilder.build.getCommandsList.asScala.map(executeCommand(_, outputDirs)).foreach {result => responseBatch.addResult(result)}
            self ! WorkComplete(JsonFormat.printToString(responseBatch.build()))
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

  private def executeCommand(command: RequestCommand, outputDirs: ListBuffer[String]): Builder = {
    def substitutePaths(path: String): String = {
      if (path.contains("{")) {
        val start = path.indexOf("{")
        val stop = path.indexOf("}", start)
        val index = path.substring(start + 1, stop).toInt
        if (index < outputDirs.size) {
          return new File(outputDirs(index), path.substring(stop + 1)).toString
        }
      }
      path
    }

    def buildCommand(outDir: String): String = {
      val executable = new File(toolsRoot, command.getExecutable).toString
      val params = command.getInput.getParamsList.asScala.map(p => s"--${p.getName} ${substitutePaths(p.getValue)}").mkString(" ")

      val cmd = s"python $executable $params --tools_root=${toolsRoot} --outdir=$outDir"
      log.debug(cmd)
      cmd
    }

    def formatOutput(response: String, outDir: String) : Builder = {
      val responseBuilder = ResponseCommand.newBuilder()
      if (response.contains("errors: none")) responseBuilder.setStatus("ok") else responseBuilder.setStatus("failed")
      responseBuilder.setMessage(response)

      for (file <- new File(outDir).listFiles()) {
        val filesBuilder = Files.newBuilder()
        filesBuilder.setName(file.getName)
        filesBuilder.setPath(file.getPath)
        responseBuilder.addFiles(filesBuilder)
      }

      responseBuilder
    }

    val group = command.getInput.getGroup
    //runNumber == subdir name in group. Name is just  an increasing number. So, latestModified dir == dir with largest name
    var runNumber = 0

    try {
      runNumber = FileUtils.getBiggestNameDir(new File(storageRoot, group).toString).toInt + 1
    } catch {
      case _: NumberFormatException => log.debug("No numeric folders in " + new File(storageRoot, group).toString)
      case _: FileNotFoundException => log.debug("No last modified folder in " + new File(storageRoot, group).toString)
    }

    val outDir = new File(storageRoot, new File(group, runNumber.toString).toString).toString
    outputDirs += outDir
    FileUtils.createDirIfNotExists(outDir)

    val outputSb = new StringBuilder
    val errorSb = new StringBuilder

    val logger = ProcessLogger((o: String) => outputSb.append(o), (e: String) => errorSb.append(e))

    Process(buildCommand(outDir), new java.io.File(toolsRoot)) ! logger
    FileUtils.scanDir(outDir, toolsRoot, group, runNumber.toString)

    val execResult = s"Result: ${outputSb.toString}, errors: ${if (errorSb.length > 0) errorSb.toString() else "none"}"
    FileUtils.saveToFile(new File(outDir, "description.txt").toString, execResult)
    val response = formatOutput(execResult, outDir)
    log.debug("Response: " + JsonFormat.printToString(response.build()))
    response
  }
}

