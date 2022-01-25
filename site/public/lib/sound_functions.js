//////////////////////////////////
// Setup
var fs = NaN;
var f0_lo = 261.63; // C4
var f0_hi = 783.99; // G5

var AudioContext = window.AudioContext || window.webkitAudioContext;
var aud_ctx;
var source;

// Stereo
var n_channels = 2;

// IDs and PCM data for piano recordings
var note_ids_C4 = ["C4", "Eb4", "E4", "G4", "C5"];
var note_pcm_C4 = new Array(5);

var note_ids_G5 = ["G5", "Bb5", "B5", "D6", "G6"];
var note_pcm_G5 = new Array(5);

//Initialize the audio context and set sample rate
function reset() {
  aud_ctx = new AudioContext();
  fs = aud_ctx.sampleRate;
}



//////////////////////////////////
// Sound generating functions

//Create array containing sequence of frequencies from an array containing a sequence of semitone distances from tonic, assuming TET    
function seq_2_freqs(seq, f0) {
  let result = new Array(seq.length);
  for (let k = 0; k < seq.length; k++) {
    result[k] = f0 * Math.pow(2, seq[k] / 12);
  }
  return result;
}



//Create a tone scramble waveform using the piano recordings
//  0 : slow, piano, G5
//  1 : slow, pure,  G5
//  2 : fast, pure,  G5
//  3 : slow, piano, C4
//  4 : slow, pure,  C4
//  5 : fast, pure,  C4
function seq_2_wave(seq, cond) {

  // Parameters
  // Note that only conditions 2 and 5 use fast pips
  let pip_dur  = (cond % 3 == 2) ? 0.0650 : 0.0650 * 32 / 12;
  let pip_len  = Math.ceil(fs * pip_dur);
  let ramp_dur = 0.0225;
  let r = Math.ceil(ramp_dur * fs);

  // Create the ramp-damp mask
  let damp = sv_prod(1 / 2, sv_sum(1, cos_vec(sv_prod(Math.PI / r, Array.from(Array(r).keys())))));
  let ramp = (cond % 3 == 0) ? ones(damp.length) : sv_sum(1, sv_prod(-1, damp)); //No ramp on for piano
  let mask = ramp.concat(ones(pip_len - 2 * r).concat(damp));

  // Create the scramble waveform
  var waveform = new Array(2);
  waveform[0] = new Array;
  waveform[1] = new Array;

  // Piano conditions (cond == 0 or 3)
  if (cond % 3 == 0) {
    for (let k = 0; k < seq.length; k++) {
      let nid = NaN;
      switch (seq[k]) {
        case 0:
          nid = 0;
          break;
        case 3:
          nid = 1;
          break;
        case 4:
          nid = 2;
          break;
        case 7:
          nid = 3;
          break;
        case 12:
          nid = 4;
          break;
      }
      
      if (cond == 0) {
        // TODO: the masking can probably be done outside of the function to reduce operations
        waveform[0] = waveform[0].concat(ew_prod(mask, v_i(note_pcm_G5[nid][0], consec(mask.length))));
        waveform[1] = waveform[1].concat(ew_prod(mask, v_i(note_pcm_G5[nid][1], consec(mask.length))));
      } else {
        // TODO: the masking can probably be done outside of the function to reduce operations
        waveform[0] = waveform[0].concat(ew_prod(mask, v_i(note_pcm_C4[nid][0], consec(mask.length))));
        waveform[1] = waveform[1].concat(ew_prod(mask, v_i(note_pcm_C4[nid][1], consec(mask.length))));
      }
    }
    return waveform;
  }
  
  // Pure-tone conditions (cond == 1, 2, 4 or 5)
  else {
    let freqs = seq_2_freqs(seq, cond > 3 ? f0_lo : f0_hi);
    for (let f = 0; f < freqs.length; f++) {
      if (Number.isNaN(freqs[f])) {
        waveform[0] = waveform[0].concat(sv_prod(0, ones(pip_len)));
      } else {
        waveform[0] = waveform[0].concat(ew_prod(mask, cos_vec(sv_prod(2 * Math.PI * freqs[f] / fs, Array.from(Array(pip_len).keys())))))
      }
    }
    waveform[1] = waveform[0];
    return waveform;
  }
}



// Create a sequence for the tone scramble
function gen_seq(cond, stim_type) {

  // Number of pips per frequency
  // One issue about comparing Exp 1 to Exp 2:
  // Exp 2 uses 32-pip scrambles in the fast condition, but Exp 1 used 12-pip scrambles in the fast conditions.
  let n_each = (cond % 3 == 2) ? 8 : 3;

  // Create the sequence (ordered completely at random)
  var SEQ = zeros(n_each).concat(sv_prod(2 + stim_type, ones(n_each))).concat(sv_prod(7, ones(n_each))).concat(sv_prod(12, ones(n_each)));
  SEQ = v_i(SEQ, randperm(SEQ.length));

  // Done
  return SEQ;
}



// Play scramble of given type
function play_scramble(scramble_type, condition_code) {
  //Create audio context if it hasn't been made already (some browsers will prevent the audio context from being created before user input)
  if (!aud_ctx) {
    reset();
  }

  //Create the tone scramble waveform / PCM data
  let seq = gen_seq(condition_code, scramble_type);
  let scramble = seq_2_wave(seq, condition_code);

  //Initialize the buffer
  let frame_count = scramble[0].length;
  let arr_buf = aud_ctx.createBuffer(n_channels, frame_count, aud_ctx.sampleRate);

  //Fill the buffer with the tone scramble waveform
  for (let channel = 0; channel < n_channels; channel++) {
    //This gives us the actual array that contains the data
    let now_buffering = arr_buf.getChannelData(channel);
    for (let k = 0; k < frame_count; k++) {
      now_buffering[k] = scramble[channel][k];
    }
  }

  //Get an AudioBufferSourceNode; this is the AudioNode to use when we want to play an AudioBuffer
  source = aud_ctx.createBufferSource();
  source.buffer = arr_buf; //Set the buffer in the AudioBufferSourceNode
  source.connect(aud_ctx.destination); //Connect the AudioBufferSourceNode to the destination so we can hear the sound
  source.start(); //Start the source playing

  //Print a message to the console when the sound is done playing
  /*source.onended = () => {
      console.log('Sound finished.');
  }*/

  //Playback initiated, return an object with all the scramble params
  var trial_obj = new Object();
  trial_obj.seq = seq;
  return trial_obj;
}



//End playback of tone scramble if subject responds early
function stop_scramble() {
  source.stop();
}



// Play calibration noise ()
function play_calib_noise() {
  // Create audio context if it hasn't been made already (some browsers will prevent the audio context from being created before user input)
  if (!aud_ctx) {
    reset();
  }

  // Create the Gaussian white noise waveform / PCM data
  // Make it three seconds long (arbitrary)
  let wn = white_noise(fs * 3)

  // Initialize the buffer
  let frame_count = wn.length;
  let arr_buf = aud_ctx.createBuffer(n_channels, frame_count, aud_ctx.sampleRate);

  // Fill the buffer with the white noise waveform
  for (let channel = 0; channel < n_channels; channel++) {
    // This gives us the actual array that contains the data
    let now_buffering = arr_buf.getChannelData(channel);
    for (let k = 0; k < frame_count; k++) {
      now_buffering[k] = wn[k];
    }
  }

  // Get an AudioBufferSourceNode; this is the AudioNode to use when we want to play an AudioBuffer
  source = aud_ctx.createBufferSource();
  source.buffer = arr_buf; //Set the buffer in the AudioBufferSourceNode
  source.connect(aud_ctx.destination); //Connect the AudioBufferSourceNode to the destination so we can hear the sound
  source.start(); //Start the source playing
}



// Play a headphone calibration stimulus
function play_calib_noise() {
  // Create audio context if it hasn't been made already (some browsers will prevent the audio context from being created before user input)
  if (!aud_ctx) {
    reset();
  }

  // Create the Gaussian white noise waveform / PCM data
  // Make it three seconds long (arbitrary)
  let wn = white_noise(fs * 3)

  // Initialize the buffer
  let frame_count = wn.length;
  let arr_buf = aud_ctx.createBuffer(n_channels, frame_count, aud_ctx.sampleRate);

  // Fill the buffer with the white noise waveform
  for (let channel = 0; channel < n_channels; channel++) {
    // This gives us the actual array that contains the data
    let now_buffering = arr_buf.getChannelData(channel);
    for (let k = 0; k < frame_count; k++) {
      now_buffering[k] = wn[k];
    }
  }

  // Get an AudioBufferSourceNode; this is the AudioNode to use when we want to play an AudioBuffer
  source = aud_ctx.createBufferSource();
  source.buffer = arr_buf; //Set the buffer in the AudioBufferSourceNode
  source.connect(aud_ctx.destination); //Connect the AudioBufferSourceNode to the destination so we can hear the sound
  source.start(); //Start the source playing
}



// Play a headphone-check stimulus
// Matches the following specifications: http://mcdermottlab.mit.edu/papers/Woods_etal_2017_headphone_screening.pdf
function play_3afc() {
  // Create audio context if it hasn't been made already (some browsers will prevent the audio context from being created before user input)
  if (!aud_ctx) {
    reset();
  }
  
  // Coding for the three types of stimuli in a trial
  //   0 : quiet target
  //   1 : in-phase distractor
  //   2 : anti-phase distractor
  let order = randperm(3);
  
  // Create the scramble waveform
  var waveform = new Array(2);
  waveform[0] = new Array;
  waveform[1] = new Array;
  
  // Silence (500 ms ISI)
  var silence = zeros(Math.floor(0.5 * fs));
  
  // Create the ramp-damp mask
  let tone_len = 1 * fs;
  let ramp_dur = 0.0225;
  let r = Math.ceil(ramp_dur * fs);
  
  let damp = sv_prod(1 / 2, sv_sum(1, cos_vec(sv_prod(Math.PI / r, Array.from(Array(r).keys())))));
  let ramp = sv_sum(1, sv_prod(-1, damp));
  let mask = ramp.concat(ones(tone_len - 2 * r).concat(damp));
  
  // Sine wave (200 Hz, 1 second duration)
  var test_tone = ew_prod(mask, cos_vec(sv_prod(2 * Math.PI * 200 / fs, Array.from(Array(tone_len).keys()))));
  
  // Generate the waveform
  for (let k = 0; k < 3; k++) {
    
    // Add ISI between each interval
    if (k > 0) {
      waveform[0] = waveform[0].concat(silence);
      waveform[1] = waveform[1].concat(silence);
    }
    
    // Add the types of tones according to the specified order
    switch (order[k]) {
      case 0:
        // Multiplying by 0.5 achieves an approx -6 dB change in intensity
        waveform[0] = waveform[0].concat(sv_prod(0.5, test_tone));
        waveform[1] = waveform[1].concat(sv_prod(0.5, test_tone));
        break;
      case 1:
        waveform[0] = waveform[0].concat(test_tone);
        waveform[1] = waveform[1].concat(test_tone);
        break;
      case 2:
        waveform[0] = waveform[0].concat(test_tone);
        waveform[1] = waveform[1].concat(sv_prod(-1, test_tone));
        break;
    }
  }
  
  // Initialize the buffer
  let frame_count = waveform[0].length;
  let arr_buf = aud_ctx.createBuffer(n_channels, frame_count, aud_ctx.sampleRate);

  // Fill the buffer with the white noise waveform
  for (let channel = 0; channel < n_channels; channel++) {
    // This gives us the actual array that contains the data
    let now_buffering = arr_buf.getChannelData(channel);
    for (let k = 0; k < frame_count; k++) {
      now_buffering[k] = waveform[channel][k];
    }
  }

  // Get an AudioBufferSourceNode; this is the AudioNode to use when we want to play an AudioBuffer
  source = aud_ctx.createBufferSource();
  source.buffer = arr_buf; //Set the buffer in the AudioBufferSourceNode
  source.connect(aud_ctx.destination); //Connect the AudioBufferSourceNode to the destination so we can hear the sound
  source.start(); //Start the source playing
  
  // Return the order of stimulus types
  return order;
}