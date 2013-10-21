import com.googlecode.protobuf.format.JsonFormat
import org.scalatest.matchers.MustMatchers
import org.scalatest.FunSpec
import protocol.Command.{ResponseCommand, RequestCommand}

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 26.09.13
 * Time: 9:23
 * To change this template use File | Settings | File Templates.
 */
class CommandProtocolTest extends FunSpec with MustMatchers {
  describe("Protocol Command") {
    it("""should parse "model list" task""") {
      val builder = RequestCommand.newBuilder()
      val jsonFormat = "{\"task\":\"model list\",\"input\":{\"group\":\"regions\"}}"
      JsonFormat.merge(jsonFormat, builder)
      val requestCommand = builder.build()
      requestCommand must have (
        'task ("model list")
      )
      requestCommand.getInput() must have (
        'group ("regions")
      )
      jsonFormat.replace(" ", "") must be (JsonFormat.printToString(requestCommand).replace(" ", ""))
    }

    it("""should parse "find pattern" and "generate model" task""") {
      val builder = RequestCommand.newBuilder()
      val jsonFormat =
        "{\"task\":\"find patterns\"," +
          "\"input\":{" +
            "\"files\":[\"file1.fasta\",\"file2.fasta\"]," +
            "\"params\":{" +
              "\"algo\" : \"random forest\"," +
              "\"algoParams\" : \"-l 10 -S 0\","+
              "\"mlWindowsize\":\"10\"," +
              "\"avgWidowsize\":\"8\"," +
              "\"mergeThreshold\":\"5\"," +
              "\"modelName\":\"/path/to/model\"" +
            "}," +
            "\"comment\":\"I am cool!\"," +
            "\"group\":\"regions\"" +
          "}," +
          "\"output\":{\"outdir\":\"/some/dir\"}" +
        "}"
      JsonFormat.merge(jsonFormat, builder)
      val requestCommand = builder.build()
      jsonFormat.replace(" ", "") must be (JsonFormat.printToString(requestCommand).replace(" ", ""))

      requestCommand must have (
        'task ("find patterns")
      )

      val input = requestCommand.getInput
      input.getFiles(0) must equal ("file1.fasta")
      input.getFiles(1) must equal ("file2.fasta")
      input must have (
        'comment ("I am cool!"),
        'group ("regions")
      )

      input.getParams must have (
        'mlWindowsize ("10"),
        'avgWidowsize ("8"),
        'mergeThreshold ("5"),
        'modelName ("/path/to/model"),
        'algo ("random forest"),
        'algoParams ("-l 10 -S 0")
      )

      requestCommand.getOutput must have (
        'outdir ("/some/dir")
      )
    }

    it("""should parse responce message""") {
      val builder = ResponseCommand.newBuilder()
      val jsonFormat =
        "{\"status\":\"ok\"," +
          "\"data\":[" +
            "{\"fullpath\":\"path\",\"description\":\"V.fasta V.kabat ...\"}," +
            "{\"fullpath\":\"path2\",\"description\" : \"V.fasta V.kabat ...\"}" +
          "]" +
        "}"
      JsonFormat.merge(jsonFormat, builder)
      val responseCommand = builder.build()
      jsonFormat.replace(" ", "") must be (JsonFormat.printToString(responseCommand).replace(" ", ""))

      responseCommand must have (
        'status ("ok")
      )

      val paths = responseCommand.getDataList

      paths.get(0) must have (
        'fullpath ("path"),
        'description ("V.fasta V.kabat ...")
      )
      paths.get(1) must have (
        'fullpath ("path2"),
        'description ("V.fasta V.kabat ...")
      )
    }
  }
}
