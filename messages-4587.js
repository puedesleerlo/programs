(Generators,socket) => {return(
    Generators.queue(notify => {
        notify(0)
      const handle = event => {
          let data = JSON.parse(event.data)
        if(data.type =="sensor") {
          return notify(data.content.variables[0].value)
        }
        }
      socket.addEventListener(`message`, handle);
      return () => socket.removeEventListener(`message`, handle);
    })
)};