import com.googlecode.protobuf.format.JsonFormat
import org.scalatest.matchers.MustMatchers
import org.scalatest.FunSpec
import protocol.Command.{ResponseCommand, BatchCommand}

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 26.09.13
 * Time: 9:23
 * To change this template use File | Settings | File Templates.
 */
class CommandProtocolTest extends FunSpec with MustMatchers {
  describe("Protocol Command") {
    it("""should parse an empty batch""") {
      val batchCommandBuilder = BatchCommand.newBuilder()
      val jsonFormat = """{"commands":[]}"""
      JsonFormat.merge(jsonFormat, batchCommandBuilder)
      val batchCommand = batchCommandBuilder.build()
      val commands = batchCommand.getCommandsList
      commands.isEmpty must be (true)
    }

    it("""should parse one simple command""") {
      val batchCommandBuilder = BatchCommand.newBuilder()
      val jsonFormat =
        """{
             "commands":[
               {
                 "executable": "ls",
                 "input": {
                     "params": [
                       {"name": "n1", "value": "v1"}, {"name": "n2", "value": "v2"}
                     ],
                     "comment": "ig is cool!",
                     "group": "all stars"
                 }
               }
             ]
          }"""
      JsonFormat.merge(jsonFormat, batchCommandBuilder)
      val batchCommand = batchCommandBuilder.build()
      val commands = batchCommand.getCommandsList
      commands.size() must be (1)
      val command = commands.get(0)
      command must have (
        'executable ("ls")
      )

      val input = command.getInput
      input must have (
        'comment ("ig is cool!"),
        'group ("all stars")
      )

      val paramsList = input.getParamsList
      paramsList.size must be (2)
      paramsList.get(0).getName must be ("n1")
      paramsList.get(0).getValue must be ("v1")
      paramsList.get(1).getName must be ("n2")
      paramsList.get(1).getValue must be ("v2")

      jsonFormat.replace(" ", "").replace("\n", "") must be (JsonFormat.printToString(batchCommand).replace(" ", ""))
    }

    it("""should parse complete batch for train and predict""") {
      val batchCommandBuilder = BatchCommand.newBuilder()
      val jsonFormat =
        """{
             "commands":[
               {
                 "executable": "train.py",
                 "input": {
                     "params": [
                       {"name": "fasta", "value": "file1.fasta"},
                       {"name": "kabat", "value": "file1.kabat"},
                       {"name": "outdir", "value": "/tmp"},
                       {"name": "model_name", "value": "model.mdl"},
                       {"name": "ml_window_size", "value": "5"}
                     ],
                     "comment": "ig is really cool!",
                     "group": "all stars"
                 }
               },
               {
                  "executable": "predict.py",
                  "input": {
                      "params": [
                        {"name": "fasta", "value": "file2.fasta"},
                        {"name": "outdir", "value": "/tmp"},
                        {"name": "model_path", "value": "/tmp/model.mdl"},
                        {"name": "merge_threshold", "value": "1"},
                        {"name": "avg_window_size", "value": "10"},
                        {"name": "ml_window_size", "value": "5"}
                      ],
                      "comment": "ig is cool!",
                      "group": "all stars"
                  }
                }
             ]
          }"""

      JsonFormat.merge(jsonFormat, batchCommandBuilder)
      val batchCommand = batchCommandBuilder.build()
      val commands = batchCommand.getCommandsList
      commands.size() must be (2)
      val trainCommand = commands.get(0)
      val predictCommand = commands.get(1)

      {
        trainCommand must have (
          'executable ("train.py")
        )

        val input = trainCommand.getInput
        input must have (
          'comment ("ig is really cool!"),
          'group ("all stars")
        )

        val paramsList = input.getParamsList
        paramsList.size must be (5)
        paramsList.get(0).getName must be ("fasta")
        paramsList.get(0).getValue must be ("file1.fasta")
        paramsList.get(1).getName must be ("kabat")
        paramsList.get(1).getValue must be ("file1.kabat")
        paramsList.get(2).getName must be ("outdir")
        paramsList.get(2).getValue must be ("/tmp")
        paramsList.get(3).getName must be ("model_name")
        paramsList.get(3).getValue must be ("model.mdl")
        paramsList.get(4).getName must be ("ml_window_size")
        paramsList.get(4).getValue must be ("5")
      }

      {
        predictCommand must have (
          'executable ("predict.py")
        )

        val input = predictCommand.getInput
        input must have (
          'comment ("ig is cool!"),
          'group ("all stars")
        )

        val paramsList = input.getParamsList
        paramsList.size must be (6)
        paramsList.get(0).getName must be ("fasta")
        paramsList.get(0).getValue must be ("file2.fasta")
        paramsList.get(1).getName must be ("outdir")
        paramsList.get(1).getValue must be ("/tmp")
        paramsList.get(2).getName must be ("model_path")
        paramsList.get(2).getValue must be ("/tmp/model.mdl")
        paramsList.get(3).getName must be ("merge_threshold")
        paramsList.get(3).getValue must be ("1")
        paramsList.get(4).getName must be ("avg_window_size")
        paramsList.get(4).getValue must be ("10")
        paramsList.get(5).getName must be ("ml_window_size")
        paramsList.get(5).getValue must be ("5")
      }
      jsonFormat.replace(" ", "").replace("\n", "") must be (JsonFormat.printToString(batchCommand).replace(" ", ""))
    }

    it("""should parse response message""") {
      val builder = ResponseCommand.newBuilder()
      val jsonFormat =
        """{"status": "ok", "message": "test"}"""
      JsonFormat.merge(jsonFormat, builder)
      val responseCommand = builder.build()
      jsonFormat.replace(" ", "").replace("\n", "") must be (JsonFormat.printToString(responseCommand).replace(" ", ""))

      responseCommand must have (
        'status ("ok"),
        'message ("test")
      )
    }
  }
}
