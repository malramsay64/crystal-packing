import { setup_state, setup_opt } from "packing_wasm"

const playPauseButton = document.getElementById("play-pause");
const image = document.getElementById("output-figure");
const stepSize = document.getElementById("step-size");
const reset = document.getElementById("reset-state");
const totalSteps = document.getElementById("total-steps");
const temperature = document.getElementById("temperature");
const coolButton = document.getElementById("cool");

let state = setup_state();
totalSteps.value = 0n;
stepSize.value = 0.1
temperature.value = Math.log10(0.2);
const optimiser = setup_opt(0.2, stepSize.value, 200n);

let cool = false;

let animationId = null;

image.innerHTML = state.as_svg();

function renderLoop () {

    optimiser.step_size = stepSize.value;
    optimiser.kt = Math.pow(10, Number(temperature.value));

    // Update state of the game
    state = optimiser.optimise_state(state);

    image.innerHTML = state.as_svg();
    stepSize.value = optimiser.step_size;
    totalSteps.value = BigInt(totalSteps.value) + optimiser.steps;

    if (cool) {
        optimiser.kt *= 0.9;
    }
    console.log(optimiser.kt);

    temperature.value = Math.log10(optimiser.kt);

    // await new Promise(r => setTimeout(r, 200));

    // Request the next frame
    animationId = requestAnimationFrame(renderLoop);
};

const isPaused = () => {
    return animationId === null;
}

const play = () => {
  playPauseButton.textContent = "⏸";
  renderLoop();
};

const pause = () => {
  playPauseButton.textContent = "▶";
  cancelAnimationFrame(animationId);
  animationId = null;
};

playPauseButton.addEventListener("click", event => {
  if (isPaused()) {
    play();
  } else {
    pause();
  }
});

coolButton.addEventListener("click", event => {
  if (!cool) {
      cool = true
  }
});

play();
