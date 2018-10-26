// var chai   = require("chai");
// var task_runner = require("../services/tasks/task_runner");
//
// describe("Scheduled Task test", function() {
//
//     before(function(done){
//         setTimeout(function(){
//             foo = true;
//             done();
//         }, 1000);
//     });
//
//     it("task test", function (done) {
//
//         var job = task_runner.schedule_task('* * * * * *', function () {
//             console.log('EXE task');
//         }, function(){
//             console.log('Complete');
//         })
//
//         job.start();
//         setTimeout(function () {
//             job.stop();
//             console.log('job stopped');
//             done();
//         }, 10000);
//     });
// });
