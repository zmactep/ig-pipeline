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
        case "find patterns" | "1" => self ! WorkComplete(findPattern(requestCommand.getInput, requestCommand.getOutput))
        case "generate model" | "2" => self ! WorkComplete(generateModel(requestCommand.getInput, requestCommand.getOutput))
        case "model list" | "3" => self ! WorkComplete(listModels(requestCommand.getInput))
        case "test" => self ! WorkComplete("{\"status\": \"ok\"}")
      }
    } pipeTo self
  }

  private def generateModel(input: RequestCommand.Input, output: RequestCommand.Output) : String = {

    def buildCommand(input: RequestCommand.Input, output: RequestCommand.Output) : String = {
      val params = Map("--ml_window_size" -> input.getParams.getMlWindowsize,
                       "--fasta" -> (if (input.getFiles(0).endsWith("fasta")) input.getFiles(0) else input.getFiles(1)),
                       "--kabat" -> (if (input.getFiles(0).endsWith("kabat")) input.getFiles(0) else input.getFiles(1)),
                       "--outdir" -> new File(storageRoot, output.getOutdir).toString,
                       "--tools_root" -> toolsRoot,
                       "--model_name" -> input.getParams.getModelName)

      val cmd = "python " + toolsRoot + "./ig-snooper/train.py " + params.map(p => p._1 + '=' + p._2).mkString(" ")
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
      if (response.contains("Your model is in ")) responseBuilder.setStatus("ok") else responseBuilder.setStatus("failed")
      val pathBuilder = PathAndDescription.newBuilder()
      pathBuilder.setDescription(response)
      pathBuilder.setFullpath(new File(new File(storageRoot, output.getOutdir).toString, input.getParams.getModelName).toString)
      responseBuilder.addData(pathBuilder.build())
      JsonFormat.printToString(responseBuilder.build())
    }

    FileUtils.createDirIfNotExists(storageRoot + output.getOutdir)

    val execResult = Process(buildCommand(input, output), new java.io.File(context.system.settings.config.getString("ig-backend.tools_root"))).!!
    saveDescription(new File(new File(storageRoot, output.getOutdir).toString, "description.txt").toString, execResult)

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
        pathBuilder.setFullpath(new File(storageRoot, model.replace("description.txt", "")).toString)
        pathBuilder.setDescription(io.Source.fromFile(new File(storageRoot, model).toString).mkString)
        responseBuilder.addData(pathBuilder.build())
      }
    }
    JsonFormat.printToString(responseBuilder.build())
  }

  private def findPattern(input: RequestCommand.Input, output: RequestCommand.Output) : String = {
    def buildCommand(input: RequestCommand.Input, output: RequestCommand.Output) : String = {
      val inputKabat = if (input.getFiles(0).endsWith("kabat")) input.getFiles(0) else input.getFiles(1)

      val params = scala.collection.mutable.Map("--ml_window_size" -> input.getParams.getMlWindowsize,
        "--fasta" -> (if (input.getFiles(0).endsWith("fasta")) input.getFiles(0) else input.getFiles(1)),
        "--outdir" -> new File(storageRoot, output.getOutdir).toString,
        "--tools_root" -> toolsRoot,
        "--model_path" -> new File(storageRoot, input.getParams.getModelPath).toString,
        "--ml_window_size" -> input.getParams.getMlWindowsize,
        "--merge_threshold" -> input.getParams.getMergeThreshold,
        "--avg_window_size" -> input.getParams.getAvgWidowsize)

      if (new File(inputKabat).exists()) params += ("--kabat" -> inputKabat)

      val cmd = "python " + toolsRoot + "./ig-snooper/predict.py " + params.map(p => p._1 + '=' + p._2).mkString(" ")
      log.debug(cmd)
      cmd
    }

    def formatOutput(response: String) : String = {
      val responseBuilder = ResponseCommand.newBuilder()
      if (response.contains("Done in ")) responseBuilder.setStatus("ok") else responseBuilder.setStatus("failed")
      val pathBuilder = PathAndDescription.newBuilder()
      pathBuilder.setDescription(response)
      pathBuilder.setFullpath(new File(new File(storageRoot, output.getOutdir).toString, "results_pic.txt").toString)
      responseBuilder.addData(pathBuilder.build())
      JsonFormat.printToString(responseBuilder.build())
    }

    FileUtils.createDirIfNotExists(storageRoot + output.getOutdir)
    val execResult = Process(buildCommand(input, output), new java.io.File(context.system.settings.config.getString("ig-backend.tools_root"))).!!

    val response = formatOutput(execResult)
    log.debug("Response: " + response)
    response
  }
}

