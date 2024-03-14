version 1.0

task max {
    input {
        Array[Int] integers
    }

    command <<<
        max=~{integers[0]}
        for i in ~{sep=" " integers}; do
            if [[ "$i" -gt "$max" ]]; then
                max=$i
            fi
        done
        echo $max
    >>>

    output {
        Int value = read_int(stdout())
    }

    runtime {
        docker: "tpesout/hpp_base:latest"
    }
}

task sum {
    input {
        Array[Int] integers
    }

    command <<<
        echo $((~{sep="+" integers}))
    >>>

    output {
        Int value = read_int(stdout())
    }

    runtime {
        docker: "tpesout/hpp_base:latest"
        preemptible: 1
    }
}
