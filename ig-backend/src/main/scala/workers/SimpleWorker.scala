package workers

import sys.process._
import scala.concurrent.Future
import akka.actor.{ActorPath, ActorRef, ActorLogging}
import akka.pattern.pipe
import protocol.Command.{ResponseCommand, RequestCommand}
import com.googlecode.protobuf.format.JsonFormat
import java.io.{IOException, FileWriter, BufferedWriter, File}
import protocol.Command.ResponseCommand.PathAndDescription
import utils.FileUtils


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
  val globalWorkdirRoot = FileUtils.fixPath(context.system.settings.config.getString("ig-backend.working_dir_root"))
  FileUtils.createDirIfNotExists(globalWorkdirRoot)

  val workDirRoot = globalWorkdirRoot + self.path.name.replace("$","") + "/"
  FileUtils.createDirIfNotExists(workDirRoot)
  FileUtils.createDirIfNotExists(storageRoot)

  def doWork(workSender: ActorRef, msg: Any): Unit = {
    Future {
      val builder = RequestCommand.newBuilder()
      try {
        msg match {
          case (m, _) => JsonFormat.merge(m.toString, builder)
        }
      } catch {
        case e: Exception => {
          log.debug(e.toString)
          self ! akka.actor.Status.Failure(e)
        }
      }
      val requestCommand = builder.build()
      requestCommand.getTask match {
        case "generate model" => self ! WorkComplete(generateModel(requestCommand.getInput, requestCommand.getOutput))
        case "model list" => self ! WorkComplete(listModels(requestCommand.getInput))
        case "find patterns" => self ! WorkComplete(findPattern(requestCommand.getInput, requestCommand.getOutput))
        //TODO handle this enum from ig-frontend
        case "1" => self ! WorkComplete(findPattern(requestCommand.getInput, requestCommand.getOutput))
        case "2" => self ! WorkComplete(generateModel(requestCommand.getInput, requestCommand.getOutput))
        case "3" => self ! WorkComplete(listModels(requestCommand.getInput))
      }
    } pipeTo self
  }

  private def generateModel(input: RequestCommand.Input, output: RequestCommand.Output) : String = {

    def buildCommand(input: RequestCommand.Input, output: RequestCommand.Output) : String = {
      val inputFasta = if (input.getFiles(0).endsWith("fasta")) input.getFiles(0) else input.getFiles(1)
      val inputKabat = if (input.getFiles(0).endsWith("kabat")) input.getFiles(0) else input.getFiles(1)

      val cmd = toolsRoot + "./ig-snooper/train_model.sh " + inputFasta + " " + inputKabat + " " + input.getParams.getMlWindowsize.toInt + " " + workDirRoot + " " + toolsRoot
      log.debug(cmd)
      cmd
    }

    def saveDescription(path: String, data: String) = {
      try {
        val out = new BufferedWriter(new FileWriter(new File(path)))
        out.write(data)
        out.close()
      } catch {
        case e: IOException =>
      }
    }

    def formatOutput(response: String) : String = {
      val responseBuilder = ResponseCommand.newBuilder()
      if (response.contains("Done. All files in")) responseBuilder.setStatus("ok") else responseBuilder.setStatus("failed")
      val pathBuilder = PathAndDescription.newBuilder()
      pathBuilder.setDescription(response)
      pathBuilder.setFullpath(storageRoot + output.getOutdir)
      responseBuilder.addData(pathBuilder.build())
      JsonFormat.printToString(responseBuilder.build())
    }

    FileUtils.createDirIfNotExists(workDirRoot)
    FileUtils.createDirIfNotExists(storageRoot + output.getOutdir)

    val execResult = Process(buildCommand(input, output), new java.io.File(context.system.settings.config.getString("ig-backend.tools_root"))).!!
    Seq("cp", "-r", workDirRoot, storageRoot + output.getOutdir).!
    saveDescription(storageRoot + output.getOutdir + "description.txt", execResult)

    val response = formatOutput(execResult)
    log.debug("Response: " + response)
    response
  }

  private def listModels(input: RequestCommand.Input) : String = {
    val models = Process("find . -name description.txt", new java.io.File(storageRoot)).!!
    val responseBuilder = ResponseCommand.newBuilder()
    responseBuilder.setStatus("ok")

    if (! models.isEmpty) {
      for (model <- models.split("\n")) {
        val pathBuilder = PathAndDescription.newBuilder()
        pathBuilder.setFullpath(storageRoot + model.replace("description.txt", ""))
        pathBuilder.setDescription(io.Source.fromFile(storageRoot + model).mkString)
        responseBuilder.addData(pathBuilder.build())
      }
    }
    JsonFormat.printToString(responseBuilder.build())
  }

  private def findPattern(input: RequestCommand.Input, output: RequestCommand.Output) : String = {
    def buildCommand(input: RequestCommand.Input, output: RequestCommand.Output) : String = {
      val inputFasta = if (input.getFiles(0).endsWith("fasta")) input.getFiles(0) else input.getFiles(1)
      val inputKabat = if (input.getFiles(0).endsWith("kabat")) input.getFiles(0) else input.getFiles(1)

      val params = input.getParams
      val cmd = toolsRoot + "./ig-snooper/predict.sh " + inputKabat + " " + inputFasta + " " + params.getMlWindowsize.toInt + " " +
        params.getAvgWidowsize.toInt + " " + params.getMergeThreshold.toInt + " " + workDirRoot + " " + toolsRoot + " " +
        storageRoot + params.getModelPath
      log.debug(cmd)
      cmd
    }

    def formatOutput(response: String) : String = {
      val responseBuilder = ResponseCommand.newBuilder()
      if (response.contains("Done. Result is in ")) responseBuilder.setStatus("ok") else responseBuilder.setStatus("failed")
      val pathBuilder = PathAndDescription.newBuilder()
      pathBuilder.setDescription(response)
      pathBuilder.setFullpath(workDirRoot)
      responseBuilder.addData(pathBuilder.build())
      JsonFormat.printToString(responseBuilder.build())
    }

    FileUtils.createDirIfNotExists(storageRoot + output.getOutdir)
    val execResult = Process(buildCommand(input, output), new java.io.File(context.system.settings.config.getString("ig-backend.tools_root"))).!!
    Seq("cp", "-r", workDirRoot, storageRoot + output.getOutdir).!

    val response = formatOutput(execResult)
    log.debug("Response: " + response)
    response
  }
}

